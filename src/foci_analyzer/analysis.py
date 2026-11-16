"""The image analysis as plain functions on arrays (plan.md, section 5).

All arrays are Z, Y, X; a 2D image, a projection and 2D labels have Z = 1, so every step is written once.
Sizes come in µm and are converted with the `Calibration`.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

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
def segment_threshold(img: np.ndarray, cal: Calibration, min_diameter: float) -> np.ndarray:
    """Otsu threshold on the smoothed image; touching objects are split on the distance map."""
    px = cal.pixel_size
    zf = cal.z_spacing / px if img.shape[0] > 1 else 1.0
    s = max(1.0, 0.3 / px) if cal.unit == "µm" else 1.5
    sm = ndi.gaussian_filter(img.astype(np.float32), (s / zf if img.shape[0] > 1 else 0, s, s))
    mask = sm > filters.threshold_otsu(sm[:, ::2, ::2])
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
    """Labels grown by `distance` (in µm or px); a negative distance gives every pixel to its nearest object."""
    if distance == 0 or labels.max() == 0:
        return labels
    d = np.inf if distance < 0 else distance / cal.pixel_size
    if labels.shape[0] == 1:
        return segmentation.expand_labels(labels[0], distance=d)[None]
    return segmentation.expand_labels(labels, distance=d, spacing=(cal.z_spacing / cal.pixel_size, 1, 1))


# ---------------------------------------------------------------- foci: filtering
_pool = ThreadPoolExecutor(max(2, min(8, (os.cpu_count() or 2) // 2)), thread_name_prefix="filter")


def together(*jobs):
    """The results of functions that run at the same time (the filters of scipy release the interpreter)."""
    futures = [_pool.submit(j) for j in jobs[1:]]
    return [jobs[0](), *(f.result() for f in futures)]


def sigma_px(size: float, cal: Calibration, slices: int, z_elongation: float) -> tuple[float, float, float]:
    """The Gaussian sigma of a focus with this diameter (its FWHM), per axis in pixels."""
    s = size / FWHM
    sz = s * z_elongation / cal.z_spacing if slices > 1 else 0.0
    return (sz, s / cal.pixel_size, s / cal.pixel_size)


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
    wavelet: product of two neighbouring à trous wavelet planes, which suppresses noise that is in one
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
        parts = together(*(second(axis) for axis, s in enumerate(sigma) if s >= 0.3))
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


def background(img: np.ndarray, sigma: tuple[float, float, float]) -> np.ndarray:
    """The image without its foci: a grey opening (in XY) with a window three focus diameters wide."""
    k = 2 * int(np.ceil(1.5 * FWHM * sigma[1])) + 1
    sm = ndi.gaussian_filter(img.astype(np.float32, copy=False), (0, 0.5 * sigma[1], 0.5 * sigma[1]))
    return ndi.grey_opening(sm, size=(1, k, k))


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

    Seeded watershed on S. A focus ends at half of its own peak (`half`) and never below `level`, and does
    not continue into another cell. Foci outside the size limits are removed.
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
    ok &= (cm == cells[keep][li - 1]) | (cm == 0)
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

    def centre(c):
        return np.bincount(li, weights=c * w, minlength=n + 1)[1:] / np.maximum(wsum, 1e-12)

    peak = np.zeros(n + 1)
    np.maximum.at(peak, li, r)
    px = cal.pixel_size
    return {"x": centre(x) * px, "y": centre(y) * px, "z": centre(z) * (cal.z_spacing if raw.shape[0] > 1 else 0),
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
