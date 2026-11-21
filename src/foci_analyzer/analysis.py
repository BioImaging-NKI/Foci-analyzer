"""The image analysis as plain functions on arrays.

All arrays are Z, Y, X; a 2D image, a projection and 2D labels have Z = 1, so every step is written once.
Sizes come in µm and are converted with the `Calibration`.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from itertools import pairwise

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import filters, morphology, segmentation

FWHM = 2.3548  # full width at half maximum of a Gaussian, in sigmas


@dataclass(frozen=True)
class Calibration:
    pixel_size: float  # XY; 1 when the image is not calibrated
    z_spacing: float
    unit: str  # "µm" or "px"

    def voxel(self, slices: int) -> float:
        """Area (Z = 1) or volume of one voxel."""
        return self.pixel_size**2 * (self.z_spacing if slices > 1 else 1.0)


# ---------------------------------------------------------------- projections
def edf(stack: np.ndarray) -> np.ndarray:
    """Extended depth of focus: for every pixel the slice where the image is sharpest."""
    s = stack.astype(np.float32)
    sharp = np.stack([ndi.gaussian_filter(np.abs(ndi.laplace(ndi.gaussian_filter(p, 1.0))), 3.0) for p in s])
    return np.take_along_axis(s, sharp.argmax(0)[None], 0)


def project(stack: np.ndarray, how: str) -> np.ndarray:
    """A stack as it is ("3D") or as one plane; the result is float32 with Z first."""
    if how == "3D" or stack.shape[0] == 1:
        return stack.astype(np.float32, copy=False)
    if how == "max projection":
        return stack.max(0)[None].astype(np.float32)
    if how == "sum projection":
        return stack.sum(0, dtype=np.float32)[None]
    if how == "extended depth of focus":
        return edf(stack)
    raise ValueError(f"unknown projection '{how}'")


def match_z(labels: np.ndarray, slices: int) -> np.ndarray:
    """Labels in the geometry of an image with `slices` slices (2D labels hold for all slices)."""
    if labels.shape[0] == slices:
        return labels
    if labels.shape[0] == 1:
        return np.broadcast_to(labels, (slices, *labels.shape[1:]))
    return labels.max(0)[None]


# ---------------------------------------------------------------- nuclei / cells
def segment_threshold(img: np.ndarray, cal: Calibration, min_diameter: float, factor: float = 1.0) -> np.ndarray:
    """Otsu threshold (times `factor`) on the smoothed image; touching objects are split on the distance map."""
    px = cal.pixel_size
    zf = cal.z_spacing / px if img.shape[0] > 1 else 1.0
    s = max(1.0, 0.3 / px) if cal.unit == "µm" else 1.5
    sm = ndi.gaussian_filter(img.astype(np.float32), (s / zf if img.shape[0] > 1 else 0, s, s))
    mask = sm > factor * filters.threshold_otsu(sm[:, ::2, ::2])
    mask = np.stack([ndi.binary_fill_holes(m) for m in mask])
    dist = ndi.distance_transform_edt(mask, sampling=(zf, 1, 1)).astype(np.float32)
    r = max(3.0, 0.5 * min_diameter / px)
    dist = ndi.gaussian_filter(dist, (0.25 * r / zf if img.shape[0] > 1 else 0, 0.25 * r, 0.25 * r))
    # one seed per object: maxima of the distance map that rise at least `h` above their saddle to a
    # higher maximum, so a dent in an outline does not split a nucleus
    # (found on a coarser grid: the seeds need no pixel precision and this is several times faster)
    # and in whole numbers of 1/16 pixel: sorting those is much faster than sorting floats)
    f = max(1, int(r / 4))
    coarse = np.round(np.minimum(dist[:, ::f, ::f], 4000) * 16).astype(np.uint16)
    peaks = morphology.h_maxima(coarse, round(16 * max(1.5, 0.2 * r))) > 0
    markers = np.zeros(dist.shape, np.int32)
    markers[:, ::f, ::f], _ = ndi.label(peaks, structure=np.ones((3, 3, 3)))
    return segmentation.watershed(-dist, markers, mask=mask).astype(np.int32)


def stitch_slices(labels: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    """3D objects from objects found per slice ("2.5D"): an object continues in the next slice when at least
    `threshold` of the smaller of the two lies on the other, best overlap first. (Relative to the smaller one,
    not the union: the top and bottom slices of a nucleus are small caps on a larger section.) A nucleus that
    a network splits in two in one slice (often in out-of-focus slices) stays one object: every piece that
    lies on it joins it.

    `labels` (Z, Y, X) has its own numbering in every slice; the result is numbered over the whole stack.
    """
    out = np.zeros(labels.shape, np.int32)
    if labels.shape[0] == 0:
        return out
    nxt = 1
    for z in range(labels.shape[0]):
        cur = labels[z].astype(np.int64)
        ids = np.unique(cur[cur > 0])
        lut = np.zeros(int(cur.max()) + 1, np.int32)
        if z > 0 and ids.size and out[z - 1].any():
            prev = out[z - 1].astype(np.int64)
            m = (cur > 0) & (prev > 0)
            k = int(prev.max()) + 1
            pairs, inter = np.unique(cur[m] * k + prev[m], return_counts=True)
            a, b = np.divmod(pairs, k)
            area_a = np.bincount(cur.ravel(), minlength=int(cur.max()) + 1)
            area_b = np.bincount(prev.ravel(), minlength=k)
            iou = inter / np.minimum(area_a[a], area_b[b])
            for i in np.argsort(-iou):
                if iou[i] < threshold:
                    break
                if lut[a[i]] == 0:  # each piece joins the object it lies on most
                    lut[a[i]] = b[i]
        for i in ids:
            if lut[i] == 0:
                lut[i] = nxt
                nxt += 1
            else:
                nxt = max(nxt, int(lut[i]) + 1)
        out[z] = lut[cur]
    return out


def filter_labels(labels: np.ndarray, cal: Calibration, min_diameter: float, max_diameter: float,
                  exclude_edges: bool) -> np.ndarray:
    """Remove objects that are too small, too large or on the image edge; number the rest 1..n."""
    foot = labels.max(0)
    n = int(labels.max())
    if n == 0:
        return labels.astype(np.int32)
    area = np.bincount(foot.ravel(), minlength=n + 1)
    diameter = 2 * np.sqrt(area / np.pi) * cal.pixel_size
    keep = (area > 0) & (diameter >= min_diameter) & (diameter <= max_diameter if max_diameter > 0 else True)
    if exclude_edges:
        edge = np.unique(np.concatenate([foot[0], foot[-1], foot[:, 0], foot[:, -1]]))
        keep[edge] = False
    keep[0] = False
    lut = np.zeros(n + 1, np.int32)
    lut[keep] = np.arange(1, int(keep.sum()) + 1)
    return lut[labels]


def expand(labels: np.ndarray, cal: Calibration, distance: float) -> np.ndarray:
    """Labels grown by `distance` (in µm or px), or shrunk when it is negative (also away from a touching
    neighbor); infinite gives every pixel to its nearest object."""
    if distance == 0 or labels.max() == 0:
        return labels
    d = abs(distance) / cal.pixel_size
    spacing = (cal.z_spacing / cal.pixel_size if labels.shape[0] > 1 else 1.0, 1.0, 1.0)
    if distance < 0:
        inner = (labels > 0) & ~segmentation.find_boundaries(labels, mode="inner")
        if labels.shape[0] == 1:
            depth = ndi.distance_transform_edt(np.pad(inner[0], 1))[1:-1, 1:-1][None]
        else:
            depth = ndi.distance_transform_edt(np.pad(inner, 1), sampling=spacing)[1:-1, 1:-1, 1:-1]
        return np.where(depth >= d, labels, 0).astype(labels.dtype)
    if labels.shape[0] == 1:
        return segmentation.expand_labels(labels[0], distance=d)[None]
    return segmentation.expand_labels(labels, distance=d, spacing=spacing)


def reach(general) -> float:
    """How far foci outside a cell still count for it (µm): infinite for the whole image."""
    return np.inf if general.whole_image else general.outside_distance


# ---------------------------------------------------------------- foci: filtering
_pool = ThreadPoolExecutor(max(2, min(8, (os.cpu_count() or 2) // 2)), thread_name_prefix="filter")


# strips of a large image filtered side by side (a pool of its own: these jobs are started from the filter pool)
_strips = ThreadPoolExecutor(max(2, min(16, os.cpu_count() or 2)), thread_name_prefix="strip")


def together(*jobs):
    """The results of functions that run at the same time (the filters of scipy release the interpreter)."""
    futures = [_pool.submit(j) for j in jobs[1:]]
    return [jobs[0](), *(f.result() for f in futures)]


def sigma_px(size: float, cal: Calibration, slices: int, z_elongation: float) -> tuple[float, float, float]:
    """The Gaussian sigma of a focus with this diameter (its FWHM), per axis in pixels."""
    s = size / FWHM
    sz = s * z_elongation / cal.z_spacing if slices > 1 else 0.0
    sxy = max(s / cal.pixel_size, 0.5)  # a focus is never smaller than a pixel in the image
    return (sz, sxy, sxy)


def _atrous(img: np.ndarray, levels: int, in_z: bool) -> list[np.ndarray]:
    """Wavelet planes of the à trous transform with the B3-spline kernel."""
    kernel = np.array([1, 4, 6, 4, 1], np.float32) / 16
    planes, smooth = [], img
    for j in range(levels):
        k = np.zeros(4 * 2**j + 1, np.float32)
        k[:: 2**j] = kernel
        nxt = smooth
        for axis in ((0, 1, 2) if in_z else (1, 2)):
            nxt = ndi.convolve1d(nxt, k, axis=axis, mode="mirror")
        planes.append(smooth - nxt)
        smooth = nxt
    return planes


def spot_filter(img: np.ndarray, sigma: tuple[float, float, float], method: str) -> np.ndarray:
    """An image in which foci of the given size stand out and the background is close to zero.

    DoG: difference of two Gaussians around the focus size. LoG: Laplacian of Gaussian at the focus size.
    wavelet: product of two neighboring à trous wavelet planes, which suppresses noise that is in one
    plane only (Olivo-Marin 2002).
    """
    img = img.astype(np.float32, copy=False)
    if method == "DoG":
        a, b = (tuple(0.8 * s for s in sigma), tuple(1.28 * s for s in sigma))
        ga, gb = together(lambda: ndi.gaussian_filter(img, a), lambda: ndi.gaussian_filter(img, b))
        ga -= gb
        return ga
    if method == "LoG":
        def second(axis):
            order = [0, 0, 0]
            order[axis] = 2
            return lambda: sigma[axis] ** 2 * ndi.gaussian_filter(img, sigma, order=order)

        # a thinner axis has no second derivative worth using
        parts = together(*(second(axis) for axis, s in enumerate(sigma) if s >= 0.3 or axis > 0))
        out = -parts[0]
        for part in parts[1:]:
            out -= part
        return out
    if method == "wavelet":
        j = max(1, round(float(np.log2(max(2 * sigma[1], 1.0)))))
        planes = _atrous(img, j + 1, in_z=sigma[0] >= 0.5)
        # Signed: positive only where both planes are positive (a focus), negative elsewhere. Keeping the
        # negative side gives the noise estimate the typical size of the product in places without foci.
        a, b = planes[j - 1], planes[j]
        mag = np.sqrt(np.abs(a * b))
        return np.where((a > 0) & (b > 0), mag, -mag)
    raise ValueError(f"unknown filter '{method}'")


# ---------------------------------------------------------------- multi-scale foci
SCALE_STEP = 1.35  # each filter size this much larger than the previous one


def scale_sizes(lo: float, hi: float) -> list[float]:
    """The foci sizes (µm) to filter at: from the smallest to the largest, about SCALE_STEP apart, at most 8."""
    lo, hi = min(lo, hi), max(lo, hi)
    if hi < 1.1 * lo:
        return [float(lo)]
    n = min(8, int(np.ceil(np.log(hi / lo) / np.log(SCALE_STEP))) + 1)
    return [float(x) for x in np.geomspace(lo, hi, n)]


ROUND_MIN = 0.3  # a peak is a focus only when it curves down in both directions (not an edge or a ridge):
# 0.3 keeps foci up to about 2.4 times as long as wide; an edge of a bright nucleus is far below it


def _second(G: np.ndarray, axis: int) -> np.ndarray:
    """Second difference along one axis (the edges repeated)."""
    a = np.moveaxis(G, axis, 0)
    out = np.empty_like(a)
    out[1:-1] = a[2:] - 2 * a[1:-1] + a[:-2]
    out[0], out[-1] = (a[1] - a[0], a[-2] - a[-1]) if len(a) > 1 else (0, 0)
    return np.moveaxis(out, 0, axis)


def _mixed_yx(G: np.ndarray) -> np.ndarray:
    """The mixed second difference in Y and X (zero at the edges)."""
    out = np.zeros_like(G)
    out[:, 1:-1, 1:-1] = (G[:, 2:, 2:] - G[:, 2:, :-2] - G[:, :-2, 2:] + G[:, :-2, :-2]) / 4
    return out


def _bin_xy(img: np.ndarray, b: int) -> np.ndarray:
    """The mean of b × b blocks in XY (the last rows and columns repeated to fill the blocks)."""
    z, y, x = img.shape
    yp, xp = -(-y // b) * b, -(-x // b) * b
    p = np.pad(img, ((0, 0), (0, yp - y), (0, xp - x)), mode="edge")
    return p.reshape(z, yp // b, b, xp // b, b).mean((2, 4), dtype=np.float32)


def _linear_axis(x: np.ndarray, b: int, start: int, stop: int, axis: int) -> np.ndarray:
    """Rows start..stop of an image b times larger along `axis`, by linear interpolation between the centers of
    the binned pixels."""
    c = (np.arange(start, stop) + 0.5) / b - 0.5
    i0 = np.clip(np.floor(c).astype(int), 0, x.shape[axis] - 1)
    i1 = np.minimum(i0 + 1, x.shape[axis] - 1)
    w = np.clip(c - i0, 0, 1).astype(np.float32)
    shape = [1] * x.ndim
    shape[axis] = len(c)
    w = w.reshape(shape)
    return np.take(x, i0, axis=axis) * (1 - w) + np.take(x, i1, axis=axis) * w


def _enlarge(x: np.ndarray, b: int, a: int, e: int, width: int) -> np.ndarray:
    """Rows a..e (all columns up to `width`) of the binned image x enlarged b times in Y and X."""
    return _linear_axis(_linear_axis(x, b, a, e, 1), b, 0, width, 2)


def _log_and_roundness_block(img: np.ndarray, sigma: tuple[float, float, float]):
    """The scale-normalized Laplacian of Gaussian (positive on a spot), its part in XY alone, and how round the
    spot is in XY: the smaller over the larger curvature of the Hessian (1 for a round spot, about 0 for an edge
    or a line). One Gaussian blur, then second differences: four times faster than four derivative filters."""
    G = ndi.gaussian_filter(img, sigma)
    yy = _second(G, 1)
    yy *= -sigma[1] ** 2  # minus the second derivatives: positive where the image curves down
    xx = _second(G, 2)
    xx *= -sigma[2] ** 2
    xy = _mixed_yx(G)
    xy *= -sigma[1] * sigma[2]
    planar = yy + xx
    if sigma[0] >= 0.3 and G.shape[0] > 1:
        F = _second(G, 0)
        F *= -sigma[0] ** 2
        F += planar
    else:
        F = planar.copy()
    del G
    # the two curvatures in XY: planar/2 ± d, d = sqrt(((yy - xx)/2)² + xy²); in place, to spare memory and time
    yy -= xx
    yy *= 0.5
    np.square(yy, out=yy)
    np.square(xy, out=xy)
    yy += xy
    d = np.sqrt(yy, out=yy)
    del xx, xy
    big = planar * 0.5
    big += d
    small = planar * 0.5
    small -= d
    roundness = np.full(planar.shape, -1.0, np.float32)
    np.divide(small, big, out=roundness, where=big > 0)
    return F.astype(np.float32, copy=False), planar.astype(np.float32, copy=False), roundness


def roundness_map(img: np.ndarray, sigma: tuple[float, float, float]) -> np.ndarray:
    """How round the image is around every voxel at one filter size (see `_log_and_roundness_block`), filtered
    in strips side by side."""
    img = img.astype(np.float32, copy=False)
    out = np.empty(img.shape, np.float32)
    rows, m = img.shape[1], int(np.ceil(4 * sigma[1])) + 2
    parts = 1 if img.size < 2_000_000 else max(1, min(4 * _strips._max_workers, rows // max(64, 8 * m)))
    edges = np.unique(np.linspace(0, rows, parts + 1).astype(int))

    def strip(a, e):
        lo, hi = max(0, a - m), min(rows, e + m)
        out[:, a:e] = _log_and_roundness_block(img[:, lo:hi], sigma)[2][:, a - lo:e - lo]

    list(_strips.map(lambda ab: strip(*ab), pairwise(edges)))
    return out


def multiscale_filter(img: np.ndarray, sizes: list[float], cal: Calibration, z_elongation: float,
                      cells: np.ndarray, per_cell: bool):
    """The Laplacian of Gaussian at every size, scale-normalized (Lindeberg) so that the sizes compare. The
    sizes vary in XY only: in Z the filter keeps the width of the smallest size, because in Z foci are mostly as
    large as the microscope makes them (a wider filter would average a thin focus away and favor the smallest
    size). Per voxel the size with the highest response in XY, the full response at that size, the index of the
    size and the roundness there; and the noise (robust SD inside the cells) of the full response per size, for
    the image and, with `per_cell`, per cell (index 0: the image value). Only the running maximum is kept."""
    img = img.astype(np.float32, copy=False)
    inside = cells > 0
    best = np.full(img.shape, -np.inf, np.float32)  # the XY response that chooses the size
    response = np.zeros(img.shape, np.float32)
    scale = np.zeros(img.shape, np.uint8)
    roundness = np.zeros(img.shape, np.float32)
    F = np.empty(img.shape, np.float32) if per_cell else None  # per cell: the whole response, for the noise per cell
    n = int(cells.max())
    objs = ndi.find_objects(cells) if per_cell and n else []
    noise_image, noise_cell = [], []
    sz = sigma_px(sizes[0], cal, img.shape[0], z_elongation)[0]
    binned = {}
    for k, size in enumerate(sizes):
        sigma = (sz, *sigma_px(size, cal, img.shape[0], z_elongation)[1:])
        # an image pyramid: a size whose filter is wider than 4 px is filtered on an image binned b × b in XY
        # (the filter then at least 2 px), and its result enlarged again; the result is smooth at such a size,
        # and the scale-normalized response does not depend on the binning
        b = 1
        while sigma[1] / (2 * b) >= 2 and img.shape[1] // (2 * b) >= 16 and img.shape[2] // (2 * b) >= 16:
            b *= 2
        small = None
        if b > 1:
            if b not in binned:
                binned[b] = _bin_xy(img, b)
            small = _log_and_roundness_block(binned[b], (sigma[0], sigma[1] / b, sigma[2] / b))

        samples = []

        def strip(a, e, sigma=sigma, k=k, samples=samples, b=b, small=small):  # one strip of rows a..e
            if small is None:  # with a margin wide enough for the filter
                m = int(np.ceil(4 * sigma[1])) + 2
                lo, hi = max(0, a - m), min(img.shape[1], e + m)
                f, xy, rnd = _log_and_roundness_block(img[:, lo:hi], sigma)
                f, xy, rnd = (x[:, a - lo:e - lo] for x in (f, xy, rnd))
            else:  # the binned result, enlarged
                f, xy, rnd = (_enlarge(x, b, a, e, img.shape[2]) for x in small)
            if F is not None:
                F[:, a:e] = f
            inner = inside[:, a:e]
            samples.append(f[inner][::8] if inner.any() else f.ravel()[::64])  # for the noise of the image
            up = xy > best[:, a:e]
            best[:, a:e][up] = xy[up]
            response[:, a:e][up] = f[up]
            scale[:, a:e][up] = k
            roundness[:, a:e][up] = rnd[up]

        # many strips, as many at a time as there are workers: the extra memory stays a fraction of the image
        rows = img.shape[1]
        margin = int(np.ceil(4 * sigma[1] / b)) + 2
        parts = 1 if img.size < 2_000_000 else max(1, min(4 * _strips._max_workers, rows // max(64, 8 * margin)))
        edges = np.unique(np.linspace(0, rows, parts + 1).astype(int) // b * b)
        edges[-1] = rows
        list(_strips.map(lambda ab: strip(*ab), pairwise(edges)))
        whole = robust_sd(np.concatenate(samples))
        per = np.full(n + 1, whole, np.float32)
        for i, sl in enumerate(objs, start=1):
            if sl is not None:
                per[i] = robust_sd(F[sl][cells[sl] == i])
        noise_image.append(whole)
        noise_cell.append(per)
    return response, scale, roundness, np.asarray(noise_image, np.float32), np.stack(noise_cell)


def suppress_nearby(pos: np.ndarray, height: np.ndarray, radius: np.ndarray, z_scale: float) -> np.ndarray:
    """Which candidates remain when, from the highest down, every remaining candidate removes the lower ones
    within its own radius (px; z counts `z_scale` px per slice). One focus per large focus, while small foci
    side by side stay apart. Indices, in the order given."""
    if len(pos) < 2:
        return np.arange(len(pos))
    pts = pos.astype(np.float64) * np.array([z_scale, 1.0, 1.0])
    order = np.argsort(-height, kind="stable")
    near = cKDTree(pts).query_ball_point(pts[order], radius[order])
    gone = np.zeros(len(pos), bool)
    for i, nb in zip(order, near, strict=True):
        if gone[i]:
            continue
        for j in nb:
            if j != i and height[j] <= height[i]:
                gone[j] = True
    return np.flatnonzero(~gone)


def multiscale_candidates(F, sigma, sizes, scale, roundness, n_img, n_cell, cells, cellmap, cal: Calibration,
                          z_elongation: float, scope: str):
    """Candidates of the multi-scale method: (unit, noise, floor, positions, heights, cells). The height of a
    voxel is its filter response divided by the noise at its own size. Only round peaks count (not the edge of
    a bright nucleus), and lower peaks within the radius of a higher one's size are left out (one focus per
    large focus)."""
    noise = float(np.median(n_img))
    if scope == "fixed value":
        unit, floor = 1.0, noise
    elif scope == "per cell":
        unit, floor = n_cell[scale, cells], 1.0
    else:
        unit, floor = n_img[scale], 1.0
    pos, height, cell = find_candidates(F / unit, sigma, floor, cellmap, True)
    if len(pos):
        ok = roundness[pos[:, 0], pos[:, 1], pos[:, 2]] >= ROUND_MIN
        pos, height, cell = pos[ok], height[ok], cell[ok]
    if len(pos):
        radius = 0.5 * np.asarray(sizes)[scale[pos[:, 0], pos[:, 1], pos[:, 2]]] / cal.pixel_size
        z_scale = cal.z_spacing / cal.pixel_size / z_elongation if F.shape[0] > 1 else 1.0
        keep = suppress_nearby(pos, height, np.maximum(radius, 1.0), z_scale)
        pos, height, cell = pos[keep], height[keep], cell[keep]
    return unit, noise, floor, pos, height, cell


def background(img: np.ndarray, sigma: tuple[float, float, float]) -> np.ndarray:
    """The image without its foci: a gray opening (in XY) with a window three focus diameters wide."""
    k = 2 * int(np.ceil(1.5 * FWHM * sigma[1])) + 1
    sm = ndi.gaussian_filter(img.astype(np.float32, copy=False), (0, 0.5 * sigma[1], 0.5 * sigma[1]))
    return ndi.grey_opening(sm, size=(1, k, k))


def background_large(img: np.ndarray, sigma: tuple[float, float, float]) -> np.ndarray:
    """As `background`, for a large focus size (multi-scale): on an image binned in XY so that the window is
    about 16 pixels, then enlarged again. The background is smooth at that size; the full opening would take
    tens of seconds on a large stack."""
    k = 2 * int(np.ceil(1.5 * FWHM * sigma[1])) + 1
    b = 1
    while k / (2 * b) >= 16 and img.shape[1] // (2 * b) >= 16 and img.shape[2] // (2 * b) >= 16:
        b *= 2
    if b == 1:
        return background(img, sigma)
    small = background(_bin_xy(img.astype(np.float32, copy=False), b), (sigma[0], sigma[1] / b, sigma[2] / b))
    out = np.empty(img.shape, np.float32)
    rows = img.shape[1]
    edges = np.unique(np.r_[np.arange(0, rows, max(64, rows // (2 * _strips._max_workers))), rows])

    def strip(a, e):
        out[:, a:e] = _enlarge(small, b, a, e, img.shape[2])

    list(_strips.map(lambda ab: strip(*ab), pairwise(edges)))
    return out


def robust_sd(values: np.ndarray) -> float:
    """1.4826 × the median absolute deviation; the plain SD when that is zero (many equal values)."""
    if values.size == 0:
        return 1.0
    if values.size > 2_000_000:
        values = values[:: values.size // 1_000_000]
    sd = 1.4826 * float(np.median(np.abs(values - np.median(values))))
    if sd <= 0:
        sd = float(values.std())
    return sd if sd > 0 else 1.0


def noise_unit(F: np.ndarray, cellmap: np.ndarray, scope: str) -> tuple[np.ndarray | float, float]:
    """What a height of 1 in the filtered image means: (unit per voxel or one number, noise of the image)."""
    inside = cellmap > 0
    whole = robust_sd(F[inside] if inside.any() else F.ravel())
    if scope == "fixed value":
        return 1.0, whole
    if scope == "per image" or not inside.any():
        return whole, whole
    unit = np.full(int(cellmap.max()) + 1, whole, np.float32)
    for i, sl in enumerate(ndi.find_objects(cellmap), start=1):
        if sl is not None:
            unit[i] = robust_sd(F[sl][cellmap[sl] == i])
    return unit[cellmap], whole


# ---------------------------------------------------------------- foci: candidates and extents
def find_candidates(S: np.ndarray, sigma: tuple[float, float, float], floor: float, cellmap: np.ndarray,
                    inside_only: bool, limit: int = 200_000):
    """Local maxima of S above `floor`: positions (n, 3: z, y, x), heights and cell numbers, highest first."""
    size = tuple(2 * max(1, round(s)) + 1 if s > 0 else 1 for s in sigma)
    peaks = (S == ndi.maximum_filter(S, size=size)) & (S >= floor)
    if inside_only:
        peaks &= cellmap > 0
    lab, n = ndi.label(peaks)  # a plateau is one candidate
    if n == 0:
        return np.zeros((0, 3), np.int32), np.zeros(0, np.float32), np.zeros(0, np.int32)
    where = np.nonzero(peaks)
    _, first = np.unique(lab[where], return_index=True)
    pos = np.column_stack([w[first] for w in where]).astype(np.int32)
    height = S[pos[:, 0], pos[:, 1], pos[:, 2]].astype(np.float32)
    order = np.argsort(-height)[:limit]
    pos, height = pos[order], height[order]
    return pos, height, cellmap[pos[:, 0], pos[:, 1], pos[:, 2]].astype(np.int32)


def foci_extents(S: np.ndarray, pos: np.ndarray, accept: np.ndarray, cells: np.ndarray, cellmap: np.ndarray,
                 level: float, half: bool, min_size: int, max_size: int):
    """Label image of the foci grown from the accepted candidates, and which candidates remain.

    Seeded watershed on S. A focus ends at half of its own peak (`half`) and never below `level`, and stays
    inside its own cell (`cellmap`: the cells grown by the distance outside cells). Foci outside the size
    limits are removed.
    """
    keep = np.flatnonzero(accept)
    labels = np.zeros(S.shape, np.int32)
    if keep.size == 0:
        return labels, keep
    p = pos[keep]
    seed = (p[:, 0], p[:, 1], p[:, 2])
    height = np.maximum(S[seed], 1e-6)
    ids = np.arange(1, keep.size + 1, dtype=np.int32)
    markers = np.zeros(S.shape, np.int32)
    markers[seed] = ids
    lowest = max(level, 0.5 * float(height.min())) if half else level
    lab = segmentation.watershed(-S, markers, mask=(S >= lowest) | (markers > 0))
    m = lab > 0
    li = lab[m]
    ok = np.ones(li.size, bool)
    if half:
        ok &= S[m] >= 0.5 * height[li - 1]
    cm = cellmap[m]
    ok &= cm == cells[keep][li - 1]
    lab[m] = np.where(ok, li, 0)
    lab[seed] = ids  # the peak itself always belongs to the focus
    count = np.bincount(lab.ravel(), minlength=keep.size + 1)[1:]
    good = (count >= min_size) & (count <= max_size)
    lut = np.zeros(keep.size + 1, np.int32)
    lut[1:][good] = np.arange(1, int(good.sum()) + 1)
    return lut[lab], keep[good]


# ---------------------------------------------------------------- measurements
def measure_foci(labels: np.ndarray, raw: np.ndarray, bg: np.ndarray, cal: Calibration) -> dict[str, np.ndarray]:
    """Per focus (label 1..n): position, size and intensities. Integrated intensity is above the background."""
    n = int(labels.max())
    z, y, x = np.nonzero(labels)
    li = labels[z, y, x]
    r = raw[z, y, x].astype(np.float64)
    above = np.clip(r - bg[z, y, x], 0, None)
    count = np.bincount(li, minlength=n + 1)[1:]
    w = np.maximum(above, 1e-6)
    wsum = np.bincount(li, weights=w, minlength=n + 1)[1:]

    def center(c):
        return np.bincount(li, weights=c * w, minlength=n + 1)[1:] / np.maximum(wsum, 1e-12)

    peak = np.zeros(n + 1)
    np.maximum.at(peak, li, r)
    px = cal.pixel_size
    return {"x": center(x) * px, "y": center(y) * px, "z": center(z) * (cal.z_spacing if raw.shape[0] > 1 else 0),
            "voxels": count, "size": count * cal.voxel(raw.shape[0]),
            "mean": np.bincount(li, weights=r, minlength=n + 1)[1:] / np.maximum(count, 1), "peak": peak[1:],
            "integrated": np.bincount(li, weights=above, minlength=n + 1)[1:]}


def per_label(labels: np.ndarray, values: np.ndarray | None, n: int) -> np.ndarray:
    """Sum of `values` (or the voxel count) per label 0..n."""
    if values is None:
        return np.bincount(labels.ravel(), minlength=n + 1).astype(np.float64)
    return np.bincount(labels.ravel(), weights=values.ravel(), minlength=n + 1)


def pearson_per_cell(a: np.ndarray, b: np.ndarray, cellmap: np.ndarray, n: int) -> np.ndarray:
    """Pearson correlation of two channels inside every cell (index 0 = outside)."""
    m = cellmap > 0
    c, a, b = cellmap[m], a[m].astype(np.float64), b[m].astype(np.float64)

    def s(w=None):
        return np.bincount(c, weights=w, minlength=n + 1)

    k = np.maximum(s(), 1)
    ma, mb = s(a) / k, s(b) / k
    cov = s(a * b) / k - ma * mb
    var = (s(a * a) / k - ma**2) * (s(b * b) / k - mb**2)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(var > 0, cov / np.sqrt(np.maximum(var, 1e-30)), np.nan)


def colocalize(lab_a: np.ndarray, lab_b: np.ndarray, xyz_a: np.ndarray, xyz_b: np.ndarray, method: str,
               min_overlap: int, max_distance: float) -> tuple[np.ndarray, np.ndarray]:
    """Which foci of A have a partner in B and the reverse (two boolean arrays)."""
    na, nb = len(xyz_a), len(xyz_b)
    with_b, with_a = np.zeros(na, bool), np.zeros(nb, bool)
    if na == 0 or nb == 0:
        return with_b, with_a
    if method == "overlap":
        m = (lab_a > 0) & (lab_b > 0)
        pairs, count = np.unique(lab_a[m].astype(np.int64) * (nb + 1) + lab_b[m], return_counts=True)
        pairs = pairs[count >= min_overlap]
        with_b[pairs // (nb + 1) - 1] = True
        with_a[pairs % (nb + 1) - 1] = True
    else:
        ta, tb = cKDTree(xyz_a), cKDTree(xyz_b)
        with_b[:] = [len(h) > 0 for h in ta.query_ball_tree(tb, max_distance)]
        with_a[:] = [len(h) > 0 for h in tb.query_ball_tree(ta, max_distance)]
    return with_b, with_a


def chance_colocalization(xyz_a: np.ndarray, radius: np.ndarray, n_b: int, cell_voxels: np.ndarray,
                          repeats: int, rng: np.random.Generator) -> float:
    """Expected fraction of the A foci with a B focus nearby when the B foci lie at random in the cell.

    `cell_voxels` are the positions (in µm) where a focus can be; `radius` is the distance per A focus
    within which a B focus counts.
    """
    if len(xyz_a) == 0 or n_b == 0 or len(cell_voxels) == 0:
        return np.nan
    hits = 0
    for _ in range(repeats):
        b = cell_voxels[rng.integers(0, len(cell_voxels), n_b)]
        d, _ = cKDTree(b).query(xyz_a)
        hits += int((d <= radius).sum())
    return hits / (repeats * len(xyz_a))


# ---------------------------------------------------------------- estimating the foci size and threshold
_DIRECTIONS = ((0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (-1, -1), (1, -1), (-1, 1))


def measure_size(img2d: np.ndarray, ys: np.ndarray, xs: np.ndarray, pixel_size: float) -> tuple[float, int]:
    """The typical diameter (FWHM, µm) of the foci at (ys, xs) and how many were measured: per focus the median
    width at half height above the local background in 8 directions, the median over the foci, corrected for
    the light smoothing used. Left out: foci near the edge, with a flat saturated top, with a brighter neighbor
    before half height (merged foci), with fewer than 6 directions that reach half height within 2 µm, or more
    than 2.5 times wider one way than another."""
    R = int(np.clip(round(2.0 / pixel_size), 8, 40))
    sm = ndi.gaussian_filter(img2d.astype(np.float32), 1.0)
    bg = ndi.minimum_filter(sm, 2 * R + 1)
    ys, xs = np.asarray(ys, int), np.asarray(xs, int)
    H, W = sm.shape
    ok = (ys >= R) & (xs >= R) & (ys < H - R) & (xs < W - R)
    flat_top = ndi.uniform_filter((img2d >= img2d.max()).astype(np.float32), 3) > 1.5 / 9
    ok[ok] &= ~flat_top[ys[ok], xs[ok]]
    ys, xs = ys[ok], xs[ok]
    if len(ys) == 0:
        return float("nan"), 0
    peak, base = sm[ys, xs], bg[ys, xs]
    half = base + (peak - base) / 2
    k = np.arange(R)
    reach = np.full((len(ys), len(_DIRECTIONS)), np.nan)
    for j, (dy, dx) in enumerate(_DIRECTIONS):
        prof = sm[ys[:, None] + k * dy, xs[:, None] + k * dx]
        below = prof <= half[:, None]
        first = np.argmax(below, 1)
        hit = below.any(1) & (first > 0)
        brighter = ((prof > peak[:, None]) & (k < first[:, None])).any(1)
        i = np.flatnonzero(hit & ~brighter)
        prev, v = prof[i, first[i] - 1], prof[i, first[i]]
        frac = np.where(prev > v, (prev - half[i]) / np.where(prev > v, prev - v, 1), 0)
        reach[i, j] = (first[i] - 1 + frac) * np.hypot(dy, dx)
    reach = reach[np.isfinite(reach).sum(1) >= 6]
    reach = reach[np.nanmax(reach, 1) <= 2.5 * np.maximum(np.nanmin(reach, 1), 0.5)] if len(reach) else reach
    if not len(reach):
        return float("nan"), 0
    widths = 2 * np.nanmedian(reach, 1)
    s = np.median(widths) / FWHM
    return float(np.sqrt(max(s * s - 1.0, 0.25)) * FWHM * pixel_size), len(reach)


def spread_sample(cells: np.ndarray, n: int = 500, per_cell: int = 10, seed: int = 0) -> np.ndarray:
    """Indices of at most n of the foci, chosen at random with at most per_cell from one cell (when there are
    more than n; with few cells more per cell, to still reach n): many foci in a few bright cells or clusters
    then do not decide the size."""
    cells = np.asarray(cells)
    if len(cells) <= n:
        return np.arange(len(cells))
    per_cell = max(per_cell, -(-n // len(np.unique(cells))))
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(cells))
    order = order[np.argsort(cells[order], kind="stable")]
    c = cells[order]
    start = np.r_[0, np.flatnonzero(np.diff(c)) + 1]
    rank = np.arange(len(c)) - np.repeat(start, np.diff(np.r_[start, len(c)]))
    pick = order[rank < per_cell]
    return np.sort(rng.choice(pick, n, replace=False) if len(pick) > n else pick)


def estimate_threshold(heights: np.ndarray, plateau: bool = True) -> tuple[float, str]:
    """A threshold (in the units of the heights) and how it was found.

    plateau: the number of foci against the threshold (on a log scale) flattens between the many low peaks of
    the noise and the foci (as TrueSpot does); the flattest point is taken. otsu: when there is no such plateau
    (or `plateau` is False), Otsu's method on the logarithm of the heights splits them into a low (noise) and a
    high (foci) group."""
    h = np.asarray(heights, float)
    h = h[h > 0]
    if len(h) < 10:
        return float("nan"), "too few peaks"
    ts = np.geomspace(1.5, max(h.max(), 3.0), 120)
    n = np.array([(h >= t).sum() for t in ts], float)
    ok = n >= 5
    if plateau and ok.sum() >= 10:
        slope = ndi.gaussian_filter1d(np.gradient(np.log(n[ok]), np.log(ts[ok])), 3)
        start = int(np.argmax(slope > slope[0] / 2)) if slope[0] < 0 else 0
        t = float(ts[ok][start + int(np.argmax(slope[start:]))])
        if t > 2.5:  # a plateau at the very bottom is the end of the noise, not a plateau
            return t, "plateau"
    lh = np.log(h[h > 1.5])
    if len(lh) < 10 or lh.max() - lh.min() < 1e-6:
        return float("nan"), "too few peaks"
    return float(np.exp(filters.threshold_otsu(lh))), "otsu"


def estimate_multiscale(raw: np.ndarray, labels: np.ndarray, cellmap: np.ndarray, cal: Calibration, general) -> dict:
    """Size range and threshold for the multi-scale method: the filter over a wide range of sizes (from 2.5
    pixels to 60 pixels or 5 µm), the threshold from the peak heights as for one size, and the range from the
    10th to the 90th percentile of the sizes of the foci above it, the upper end one filter step wider (the
    lower end stays: below 2.5 pixels a size cannot be told apart from the pixels)."""
    lo, hi = max(0.05, 2.5 * cal.pixel_size), min(5.0, 60 * cal.pixel_size)
    sizes = scale_sizes(lo, hi)
    cells = match_z(labels, raw.shape[0])
    F, scale, rnd, n_img, n_cell = multiscale_filter(raw, sizes, cal, general.z_elongation, cells, False)
    sigma = sigma_px(sizes[0], cal, raw.shape[0], general.z_elongation)
    _, noise, _, pos, h, _ = multiscale_candidates(F, sigma, sizes, scale, rnd, n_img, n_cell, cells, cellmap, cal,
                                                   general.z_elongation, "per image")
    t, how = estimate_threshold(h)
    if not np.isfinite(t):
        raise ValueError("too few peaks to estimate a threshold")
    t = float(min(50.0, max(1.0, round(t, 1))))
    found = np.asarray(sizes)[scale[pos[:, 0], pos[:, 1], pos[:, 2]]][h >= t]
    if len(found) < 3:
        raise ValueError("too few foci to estimate the size range")
    small = max(0.05, round(float(np.percentile(found, 10)) / 0.05) * 0.05)
    large = max(small * SCALE_STEP, round(float(np.percentile(found, 90)) * SCALE_STEP / 0.05) * 0.05)
    large = round(large / 0.05) * 0.05
    return {"size": round(small, 2), "size_max": round(large, 2), "threshold": t, "how": how, "noise": noise,
            "foci": int((h >= t).sum())}


def estimate_foci(stack: np.ndarray, labels: np.ndarray, cal: Calibration, general, fc) -> dict:
    """Foci size and threshold for one foci channel (`stack`: Z, Y, X of the prepared image) in its cells."""
    raw = project(stack, general.dimension)
    cellmap = match_z(expand(labels, cal, reach(general)), raw.shape[0])
    mask = cellmap.max(0) > 0
    if not mask.any():
        raise ValueError("no cells to estimate in")
    if fc.method == "multi-scale":
        return estimate_multiscale(raw, labels, cellmap, cal, general)
    return estimate_watershed(raw, cellmap, cal, general, fc)


def estimate_watershed(raw: np.ndarray, cellmap: np.ndarray, cal: Calibration, general, fc, passes: int = 4) -> dict:
    """Size and threshold for the watershed method, found together: the image is filtered at a size, the
    threshold estimated from the peak heights, and the size measured on the foci above that threshold (at most
    500, spread over the cells); then again at the measured size, until the size stays the same. It starts at a
    small size (0.4 µm, at least 3 pixels), at which foci of every size give a peak. A plateau that leaves too
    few foci to measure (a few very bright outliers) gives way to Otsu's method."""
    img2d = raw.max(0)
    scope = "per image" if fc.threshold_scope == "fixed value" else fc.threshold_scope
    size, used = max(0.4, round(3 * cal.pixel_size / 0.05) * 0.05), 0
    for i in range(passes):
        sigma = sigma_px(size, cal, raw.shape[0], general.z_elongation)
        F = spot_filter(raw, sigma, fc.filter)
        unit, noise = noise_unit(F, cellmap, scope)
        pos, h, cells = find_candidates(F / unit, sigma, 1.0, cellmap, True)
        for plateau in (True, False):
            t, how = estimate_threshold(h, plateau)
            if not np.isfinite(t):
                raise ValueError("too few peaks to estimate a threshold")
            foci = np.flatnonzero(h >= t)
            pick = foci[spread_sample(cells[foci])]
            measured, used = measure_size(img2d, pos[pick, 1], pos[pick, 2], cal.pixel_size)
            if (used >= 8 and np.isfinite(measured)) or how != "plateau":
                break
        if i == passes - 1:
            break
        if used < 8 or not np.isfinite(measured):
            raise ValueError("too few clear foci to measure the size of")
        measured = max(0.05, round(measured / 0.05) * 0.05)
        if abs(measured - size) < 0.025:
            break
        size = measured
    t = float(min(50.0, max(1.0, round(t, 1))))
    return {"size": round(size, 2), "threshold": t, "how": how, "noise": float(noise), "foci": int((h >= t).sum()),
            "measured": used}
