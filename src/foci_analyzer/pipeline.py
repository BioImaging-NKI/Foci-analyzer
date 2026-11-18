"""The analysis of one frame of one image as a chain of cached stages (plan.md, section 5).

prepare → segment → per foci channel: filter → candidates → extents → measure → tables per cell and focus.
A stage only runs again when its own settings or an earlier stage changed, so moving the threshold slider
repeats the extents and measurements only.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field, replace
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from . import analysis as an
from . import models
from .config import Settings
from .io import ImageInfo
from .stages import StageCache


@dataclass
class ChannelResult:
    channel: int  # 1-based
    method: str
    noise: float  # robust SD of the filtered image inside the cells
    threshold: float  # in the units of `height`
    floor: float  # lowest height that is kept as a candidate
    height_unit: str  # "× noise", "filtered" or "probability"
    pos: np.ndarray  # candidates: (n, 3) z, y, x in pixels of the detection image
    height: np.ndarray
    cell: np.ndarray
    kept: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))  # candidates that are foci
    labels: np.ndarray | None = None
    table: dict = field(default_factory=dict)
    raw: np.ndarray | None = None  # the image the foci were detected in (a stack or a projection)
    filtered: tuple | None = None  # (filtered image, divisor): filtered / divisor is in the units of `height`

    def summary(self) -> dict:
        return {"channel": self.channel, "method": self.method, "noise": self.noise, "threshold": self.threshold,
                "floor": self.floor, "height_unit": self.height_unit, "candidates": len(self.height),
                "foci": len(self.kept), "filtered": self.filtered is not None}


@dataclass
class FrameResult:
    key: str
    cal: an.Calibration
    binning: int
    border: int
    labels: np.ndarray  # nuclei or cells, Z = 1 for 2D segmentation
    channels: list[ChannelResult]
    cells: pd.DataFrame
    foci: pd.DataFrame
    remarks: list[str]
    offset: tuple[int, int] = (0, 0)  # x, y of the analyzed part in the image (preview of the visible part only)
    removed: list[int] = field(default_factory=list)  # cells removed by hand: not in the tables (curation.py)


def _prepare(data: np.ndarray, info: ImageInfo, s: Settings):
    """Channels (C, Z, Y, X) after removing the border and binning, with the calibration that then holds."""
    b, n = s.image.crop_border, s.image.xy_binning
    if b and min(data.shape[-2:]) > 2 * b + 8:
        data = data[..., b:-b, b:-b]
    else:
        b = 0
    if n > 1:
        c, z, h, w = data.shape
        h, w = h // n * n, w // n * n
        data = data[..., :h, :w].reshape(c, z, h // n, n, w // n, n).mean(axis=(3, 5), dtype=np.float32)
    px = (info.pixel_size or 1.0) * n
    cal = an.Calibration(px, info.z_spacing or px, info.unit)
    return data, cal, b


def _read_labels(folder: str, name: str, shape: tuple) -> np.ndarray:
    for ext in (".tif", ".tiff", ".png"):
        for stem in (name, f"{name}_labels", f"{name}__labelmap"):
            p = Path(folder) / f"{stem}{ext}"
            if p.is_file():
                lab = np.squeeze(tifffile.imread(p) if ext != ".png" else _png(p))
                lab = lab[None] if lab.ndim == 2 else lab
                if lab.ndim != 3 or lab.shape[1:] != shape[1:]:
                    raise ValueError(f"the label image {p.name} does not have the size of the image")
                return lab.astype(np.int32)
    raise ValueError(f"no label image for '{name}' in {folder or '(no labels folder set)'}")


def _png(path: Path) -> np.ndarray:
    from skimage import io as skio
    return skio.imread(path)


def _segment(data: np.ndarray, cal: an.Calibration, s: Settings, name: str) -> np.ndarray:
    seg = s.segmentation
    ch = s.image.cell_channel or s.image.nuclei_channel
    if ch == 0:
        return np.ones((1, *data.shape[-2:]), np.int32)  # no cells: the whole image is one
    if ch > data.shape[0]:
        raise ValueError(f"the image has no channel {ch}")
    img = data[ch - 1]
    if seg.dimension == "2D" and img.shape[0] > 1:
        img = img.max(0)[None]
    if seg.method == "threshold":
        labels = an.segment_threshold(img, cal, seg.min_diameter, seg.threshold_factor)
    elif seg.method == "labels from file":
        labels = _read_labels(seg.labels_folder, name, img.shape)
    else:
        if seg.method == "stardist" and s.image.cell_channel:
            raise ValueError("StarDist finds nuclei only: choose Cellpose or InstanSeg for a cell channel, or set the cell channel to 0")
        second = data[s.image.nuclei_channel - 1] if s.image.cell_channel and s.image.nuclei_channel else None
        if second is not None and second.shape[0] != img.shape[0]:
            second = second.max(0)[None]
        labels = models.segment(img, second, cal, seg)
        if seg.method in ("stardist", "instanseg") and labels.shape[0] > 1:
            labels = an.stitch_slices(labels)  # 2D networks: every slice, then linked into 3D objects
    return an.filter_labels(labels, cal, seg.min_diameter, seg.max_diameter, seg.exclude_on_edges)


def analyze_frame(cache: StageCache, image_key: str, data: np.ndarray, info: ImageInfo, s: Settings,
                  only_candidates: bool = False) -> FrameResult:
    """Analyze one frame (`data`: C, Z, Y, X). `image_key` names the pixels (path, change time and frame)."""
    remarks: list[str] = []
    if info.slices > 1 and info.z_spacing is None and info.unit == "µm":
        remarks.append("no Z spacing in the file: taken equal to the pixel size")
    if s.segmentation.dimension == "3D" and info.slices > 1 and s.foci_general.dimension != "3D":
        # cells in 3D and foci on a projection do not fit together: a projected focus has no slice to be in
        s = replace(s, foci_general=replace(s.foci_general, dimension="3D"))
        remarks.append("the cells are segmented in 3D, so the foci are detected in 3D as well")
    steps = 2 + 3 * len(s.foci) + (0 if only_candidates else 2)
    k_prep, (data, cal, border) = cache.run("prepare", asdict(s.image) | {"k": "p"}, (image_key,),
                                            lambda: _prepare(data, info, s), "Preparing the image", (1, steps))
    seg_params = asdict(s.segmentation) | {"n": s.image.nuclei_channel, "c": s.image.cell_channel}
    what = "cells" if s.image.cell_channel else "nuclei"
    k_seg, labels = cache.run("segment", seg_params, (k_prep,), lambda: _segment(data, cal, s, info.name),
                              f"Finding the {what} ({s.segmentation.method})", (2, steps))
    n_cells = int(labels.max())
    general = s.foci_general

    for fc in s.foci:
        if fc.channel > data.shape[0]:
            raise ValueError(f"the image has no channel {fc.channel}")

    def one(n, fc):  # one foci channel; the channels are done at the same time

        def detect(fc=fc):
            raw = an.project(data[fc.channel - 1], general.dimension)
            sigma = an.sigma_px(fc.size, cal, raw.shape[0], general.z_elongation)
            F, bg = an.together(lambda: an.spot_filter(raw, sigma, fc.filter), lambda: an.background(raw, sigma))
            return raw, F, bg, sigma

        k_det, (raw, F, bg, sigma) = cache.run(
            "filter", {"ch": fc.channel, "dim": general.dimension, "size": fc.size, "filter": fc.filter,
                       "ze": general.z_elongation}, (k_prep,), detect,
            f"Filtering channel {fc.channel}", (3 + 3 * n, steps))

        def candidates(fc=fc, raw=raw, F=F, sigma=sigma):
            cellmap = an.match_z(an.expand(labels, cal, an.reach(general)), raw.shape[0])
            unit, noise = an.noise_unit(F, an.match_z(labels, raw.shape[0]), fc.threshold_scope)
            if fc.method == "spotiflow":
                pos, prob = models.spots(raw, cal, fc)
                inside = cellmap[pos[:, 0], pos[:, 1], pos[:, 2]] > 0
                pos, prob = pos[inside], prob[inside]
                return cellmap, noise, noise, 0.05, pos, prob, cellmap[pos[:, 0], pos[:, 1], pos[:, 2]]
            floor = noise if fc.threshold_scope == "fixed value" else 1.0
            S = F if fc.threshold_scope == "fixed value" else F / unit
            return (cellmap, unit, noise, floor, *an.find_candidates(S, sigma, floor, cellmap, True))

        k_cand, (cellmap, unit, noise, floor, pos, height, cell) = cache.run(
            "candidates", {"scope": fc.threshold_scope, "out": an.reach(general), "method": fc.method,
                           "model": fc.spotiflow_model if fc.method == "spotiflow" else ""}, (k_det, k_seg),
            candidates, f"Finding foci in channel {fc.channel}"
            + (" (Spotiflow)" if fc.method == "spotiflow" else ""), (4 + 3 * n, steps))

        if fc.method == "spotiflow":
            threshold, height_unit = fc.probability, "probability"
        elif fc.threshold_scope == "fixed value":
            threshold, height_unit = fc.fixed_threshold or 5 * noise, "filtered"
        else:
            threshold, height_unit = fc.sensitivity, "× noise"
        res = ChannelResult(fc.channel, fc.method, noise, threshold, floor, height_unit, pos, height, cell)
        res.raw = raw
        if fc.method != "spotiflow":  # for the viewer (Filtered): no copy, the cache holds the filtered image anyway
            res.filtered = (F, 1.0 if fc.threshold_scope == "fixed value" else unit)
        if only_candidates:
            return res, None

        def extents(fc=fc, F=F, unit=unit, pos=pos, height=height, cell=cell, cellmap=cellmap, raw=raw, bg=bg,
                    threshold=threshold):
            S = F / unit
            if fc.method == "spotiflow":
                level, half = 1.0, True
            else:
                half = fc.extent == "half maximum"
                level = (threshold / unit if fc.threshold_scope == "fixed value" else threshold) if not half else 0.0
            lab, kept = an.foci_extents(S, pos, height >= threshold, cell, cellmap, level, half, fc.min_size,
                                        fc.max_size)
            table = an.measure_foci(lab, raw, bg, cal)
            table["cell"], table["height"] = cell[kept], height[kept]
            return lab, kept, table

        k_ext, (res.labels, res.kept, res.table) = cache.run(
            "extents", {"thr": threshold, "extent": fc.extent, "min": fc.min_size, "max": fc.max_size},
            (k_cand,), extents, f"Outlining and measuring the foci of channel {fc.channel}", (5 + 3 * n, steps))
        return res, k_ext

    if len(s.foci) > 1:
        with ThreadPoolExecutor(len(s.foci), thread_name_prefix="channel") as pool:
            done = list(pool.map(one, range(len(s.foci)), s.foci))
    else:
        done = [one(n, fc) for n, fc in enumerate(s.foci)]
    channels = [d[0] for d in done]
    keys = [k_seg, *(d[1] for d in done if d[1])]

    if only_candidates:
        empty = pd.DataFrame()
        return FrameResult("", cal, s.image.xy_binning, border, labels, channels, empty, empty, remarks)

    pairs = [c.channel for c in channels] if s.colocalization.enabled else []
    _, base = cache.run("cells", {"pairs": pairs, "dim": general.dimension}, (k_seg,),
                        lambda: _cell_base(data, labels, n_cells, cal, channels, pairs),
                        "Measuring the cells", (steps - 1, steps))

    def tables():
        return _tables(data, base, labels, n_cells, cal, channels, s)

    key, (cells, foci) = cache.run("tables", asdict(s.colocalization), tuple(keys), tables,
                                  "Making the tables", (steps, steps))
    remarks = remarks + saturation_notes(foci, data.dtype)
    return FrameResult(key, cal, s.image.xy_binning, border, labels, channels, cells, foci, remarks)


def saturation_notes(foci: pd.DataFrame, dtype) -> list[str]:
    """Per foci channel of an image: a note when more than a fifth of its foci reach the highest gray value
    (once per image rather than per cell: in 8-bit images that is often almost every cell)."""
    if not np.issubdtype(dtype, np.integer) or not len(foci) or "peak" not in foci:
        return []
    top = np.iinfo(dtype).max
    out = []
    for c, g in foci.groupby("channel"):
        share = float((g["peak"] >= top).mean())
        if share > 0.2:
            out.append(f"saturated foci in channel {int(c)}: {round(100 * share)} % reach the highest gray value "
                       "(the counts are fine, their intensities are not)")
    return out


def _cell_base(data, labels, n, cal, channels: list[ChannelResult], pairs: list[int]):
    """What does not depend on the foci: position, size and intensities per cell, Pearson per channel pair."""
    foot = labels.max(0)
    area = an.per_label(foot, None, n)
    yy, xx = np.nonzero(foot)
    li = foot[yy, xx]
    cells = pd.DataFrame({"cell": np.arange(1, n + 1)})
    k = np.maximum(area[1:], 1)
    cells["x"] = np.bincount(li, weights=xx, minlength=n + 1)[1:] / k * cal.pixel_size
    cells["y"] = np.bincount(li, weights=yy, minlength=n + 1)[1:] / k * cal.pixel_size
    cells["area"] = area[1:] * cal.pixel_size**2
    vox = an.per_label(labels, None, n)[1:]
    if labels.shape[0] > 1:
        cells["volume"] = vox * cal.voxel(labels.shape[0])
    edge = np.zeros(n + 1, bool)
    if foot.min() == 0:  # not when the whole image is the one cell
        edge[np.unique(np.concatenate([foot[0], foot[-1], foot[:, 0], foot[:, -1]]))] = True
    for c in range(data.shape[0]):
        if labels.shape[0] == 1:
            plane = data[c].mean(0, dtype=np.float32)[None] if data.shape[1] > 1 else data[c]
            cells[f"mean_ch{c + 1}"] = an.per_label(labels, plane, n)[1:] / k
        else:
            cells[f"mean_ch{c + 1}"] = an.per_label(labels, data[c], n)[1:] / np.maximum(vox, 1)
    by_channel = {c.channel: c for c in channels}
    for a, b in combinations(pairs, 2):
        ra, rb = by_channel[a].raw, by_channel[b].raw
        cells[f"pearson_ch{a}_ch{b}"] = an.pearson_per_cell(ra, rb, an.match_z(labels, ra.shape[0]), n)[1:]
    return cells, edge[1:]


def _tables(data, base, labels, n, cal, channels: list[ChannelResult], s: Settings):
    """The table per cell and the table per focus of one frame."""
    cells = base[0].copy()
    flags = [["on edge"] if e else [] for e in base[1]]
    # size much off the typical cell of this frame (from 5 cells on): two nuclei found as one, or debris
    area = cells["area"].to_numpy()
    if n >= 5:
        med = np.median(area)
        for i in np.flatnonzero(area > 2 * med):
            flags[i].append("much larger than the others (merged?)")
        for i in np.flatnonzero(area < 0.4 * med):
            flags[i].append("much smaller than the others (debris?)")

    parts = []
    for res in channels:
        t, c = res.table, res.channel
        cnt = np.bincount(t["cell"], minlength=n + 1)[1:n + 1]
        cells[f"foci_ch{c}"] = cnt
        for name, col in (("intensity", "integrated"), ("size", "size")):
            cells[f"foci_ch{c}_total_{name}"] = np.bincount(t["cell"], weights=t[col], minlength=n + 1)[1:n + 1]
        cells[f"foci_ch{c}_mean_size"] = np.where(cnt > 0, cells[f"foci_ch{c}_total_size"] / np.maximum(cnt, 1), np.nan)
        cells[f"foci_ch{c}_mean_intensity"] = np.where(cnt > 0, cells[f"foci_ch{c}_total_intensity"] / np.maximum(cnt, 1), np.nan)
        if n >= 5:  # far more foci than the typical cell: robustly above the rest, and at least twice the median
            med = np.median(cnt)
            mad = 1.4826 * np.median(np.abs(cnt - med))
            for i in np.flatnonzero((cnt > 2 * med) & (cnt > med + 4 * max(mad, 1.0))):
                flags[i].append(f"many foci in channel {c}")
        df = pd.DataFrame({"cell": t["cell"], "channel": c, "focus": np.arange(1, len(t["cell"]) + 1)})
        for col in ("x", "y", "z", "size", "voxels", "mean", "peak", "integrated", "height"):
            df[col] = t[col]
        parts.append(df)

    co = s.colocalization
    if co.enabled and len(channels) > 1:
        rng = np.random.default_rng(1)
        for (i, a), (j, b) in combinations(enumerate(channels), 2):
            ta, tb = a.table, b.table
            xyz_a, xyz_b = (np.column_stack([t["x"], t["y"], t["z"]]) for t in (ta, tb))
            with_b, with_a = an.colocalize(a.labels, b.labels, xyz_a, xyz_b, co.method, co.min_overlap,
                                           co.max_distance)
            parts[i][f"with_ch{b.channel}"], parts[j][f"with_ch{a.channel}"] = with_b, with_a
            pair = f"ch{a.channel}_ch{b.channel}"
            both = np.bincount(ta["cell"], weights=with_b, minlength=n + 1)[1:n + 1]
            both_b = np.bincount(tb["cell"], weights=with_a, minlength=n + 1)[1:n + 1]
            na, nb = cells[f"foci_ch{a.channel}"].to_numpy(), cells[f"foci_ch{b.channel}"].to_numpy()
            cells[f"coloc_{pair}"] = both.astype(int)
            cells[f"fraction_ch{a.channel}_with_ch{b.channel}"] = np.where(na > 0, both / np.maximum(na, 1), np.nan)
            cells[f"fraction_ch{b.channel}_with_ch{a.channel}"] = np.where(nb > 0, both_b / np.maximum(nb, 1), np.nan)
            raw_a, raw_b = a.raw, b.raw
            cm = an.match_z(labels, raw_a.shape[0])
            # Manders: share of the foci intensity of one channel that lies inside foci of the other
            for (p, q, rp, tp) in ((a, b, raw_a, ta), (b, a, raw_b, tb)):
                m = p.labels > 0
                cell_of = np.concatenate([[0], tp["cell"]])[p.labels[m]]
                w = rp[m].astype(np.float64)
                tot = np.bincount(cell_of, weights=w, minlength=n + 1)[1:n + 1]
                ins = np.bincount(cell_of, weights=w * (q.labels[m] > 0), minlength=n + 1)[1:n + 1]
                cells[f"manders_ch{p.channel}_in_ch{q.channel}"] = np.where(tot > 0, ins / np.maximum(tot, 1e-12), np.nan)
            if co.chance_repeats > 0:
                cells[f"chance_ch{a.channel}_with_ch{b.channel}"] = _chance(cm, cal, ta, tb, xyz_a, co, n, rng)
                cells[f"chance_ch{b.channel}_with_ch{a.channel}"] = _chance(cm, cal, tb, ta, xyz_b, co, n, rng)
    cells["check"] = ["; ".join(f) for f in flags]
    foci = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return cells, foci


def _chance(cm, cal, ta, tb, xyz_a, co, n, rng) -> np.ndarray:
    from scipy import ndimage as ndi
    out = np.full(n, np.nan)
    dim3 = cm.shape[0] > 1

    def radius(size):
        return (3 * size / (4 * np.pi)) ** (1 / 3) if dim3 else np.sqrt(size / np.pi)

    rb = float(np.median(radius(tb["size"]))) if len(tb["size"]) else 0.0
    for i, sl in enumerate(ndi.find_objects(cm), start=1):
        ia, nb = np.flatnonzero(ta["cell"] == i), int((tb["cell"] == i).sum())
        if sl is None or ia.size == 0 or nb == 0:
            continue
        z, y, x = np.nonzero(cm[sl] == i)
        vox = np.column_stack([(x + sl[2].start) * cal.pixel_size, (y + sl[1].start) * cal.pixel_size,
                               (z + sl[0].start) * (cal.z_spacing if dim3 else 0)])
        r = np.full(ia.size, co.max_distance) if co.method == "centroid distance" else radius(ta["size"][ia]) + rb
        out[i - 1] = an.chance_colocalization(xyz_a[ia], r, nb, vox, co.chance_repeats, rng)
    return out
