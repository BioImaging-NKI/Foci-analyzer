"""What the analysis leaves besides the tables: the store the viewer reads back, label images, RGB pictures and
ImageJ TIFFs with overlays.

All label images of an analysis are in the grid the analysis saw: after XY binning and without the removed
border (and, for a preview of the visible part, from its offset). `to_image()` brings coordinates back to
pixels of the image.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage import measure

from . import analysis as an
from .pipeline import ChannelResult, FrameResult

NUCLEI, CELLS = (0, 128, 255), (255, 128, 0)
OTHER = [(0, 255, 0), (255, 0, 255), (255, 255, 0), (255, 255, 255)]
OUTLINE = [(255, 224, 77), (77, 224, 255), (255, 153, 51), (204, 153, 255)]  # foci channel 1, 2, … (as the viewer)
META_COLUMNS = ["condition", "image", "folder", "path", "frame"]


def safe_name(text: str) -> str:
    return re.sub(r"[^\w.-]+", "_", text).strip("_") or "image"


def default_colours(n_channels: int, nuclei: int, cell: int, foci: list[int]) -> list[tuple[int, int, int]]:
    """The channel colours of the viewer when the user did not choose any."""
    special = [nuclei - 1, cell - 1]
    order = [c - 1 for c in dict.fromkeys(foci) if c - 1 not in special]
    order += [c for c in range(n_channels) if c not in special and c not in order]
    out = []
    for c in range(n_channels):
        out.append(NUCLEI if c == nuclei - 1 else CELLS if c == cell - 1 else OTHER[order.index(c) % len(OTHER)])
    return out


# ---------------------------------------------------------------- the store the viewer reads back
def store_path(folder: str | Path, image_id: str) -> Path:
    return Path(folder) / ".viewer" / (hashlib.sha1(image_id.encode()).hexdigest()[:16] + ".npz")


def _small(a: np.ndarray) -> np.ndarray:
    return a.astype(np.uint16) if a.size and int(a.max()) < 65536 else a.astype(np.int32)


def save_store(path: Path, frames: list[FrameResult]) -> None:
    """Everything the viewer needs to show the result of a batch run again, per frame, compressed."""
    arrays, meta = {}, []
    for t, r in enumerate(frames):
        arrays[f"f{t}_cells"] = _small(r.labels)
        chans = []
        for i, ch in enumerate(r.channels):
            p = f"f{t}_c{i}_"
            arrays[p + "labels"] = _small(ch.labels)
            arrays[p + "pos"], arrays[p + "height"] = ch.pos.astype(np.int32), ch.height.astype(np.float32)
            arrays[p + "cell"], arrays[p + "kept"] = ch.cell.astype(np.int32), np.asarray(ch.kept, np.int64)
            chans.append({k: getattr(ch, k) for k in ("channel", "method", "noise", "threshold", "floor", "height_unit")})
        meta.append({"binning": r.binning, "border": r.border, "remarks": r.remarks, "channels": chans,
                     "cal": [r.cal.pixel_size, r.cal.z_spacing, r.cal.unit]})
    arrays["meta"] = np.frombuffer(json.dumps(meta).encode(), np.uint8)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, **arrays)
    tmp.replace(path)


def load_store(path: Path, t: int, cells: pd.DataFrame, foci: pd.DataFrame) -> FrameResult | None:
    """Frame t from the store, with its rows of the result tables; None when it is not there."""
    if not path.is_file():
        return None
    with np.load(path) as f:
        meta = json.loads(bytes(f["meta"]).decode())
        if t >= len(meta):
            return None
        m = meta[t]
        channels = []
        for i, c in enumerate(m["channels"]):
            p = f"f{t}_c{i}_"
            ch = ChannelResult(c["channel"], c["method"], c["noise"], c["threshold"], c["floor"], c["height_unit"],
                               f[p + "pos"].reshape(-1, 3), f[p + "height"], f[p + "cell"])
            ch.kept, ch.labels = f[p + "kept"], f[p + "labels"].astype(np.int32)
            channels.append(ch)
        labels = f[f"f{t}_cells"].astype(np.int32)
    cal = an.Calibration(*m["cal"])
    drop = [c for c in META_COLUMNS if c in cells]
    cells = cells.drop(columns=drop).reset_index(drop=True)
    foci = foci.drop(columns=[c for c in META_COLUMNS if c in foci]).reset_index(drop=True)
    key = "stored|" + hashlib.sha1(f"{path}|{path.stat().st_mtime}|{t}".encode()).hexdigest()[:16]
    return FrameResult(key, cal, m["binning"], m["border"], labels, channels, cells, foci, m["remarks"])


# ---------------------------------------------------------------- label images
def save_labels(folder: Path, base: str, frames: list[FrameResult], tracks: dict | None = None) -> None:
    """The label images of one image as TIFF (T, Z, Y, X): the cells and the foci of every foci channel.

    They are in the grid of the analysis (after binning, without the removed border; the pixel size in the
    file says so). In a tracked time-lapse the cells are numbered by their track, so one cell keeps its number
    over time (`tracks`: per frame {cell: track}); otherwise numbers are per frame, as in the tables.
    """
    folder.mkdir(parents=True, exist_ok=True)
    cal = frames[0].cal
    cells = []
    for t, r in enumerate(frames):
        lab = r.labels
        if tracks and tracks.get(t):
            lut = np.zeros(int(lab.max()) + 1, np.int64)
            for cell, track in tracks[t].items():
                if cell < len(lut):
                    lut[cell] = track
            lab = lut[lab]
        cells.append(lab)
    stacks = {"cells": cells}
    for i, ch in enumerate(frames[0].channels):
        stacks[f"foci_ch{ch.channel}"] = [f.channels[i].labels for f in frames]
    for name, planes in stacks.items():
        depth = max(p.shape[0] for p in planes)
        stack = np.stack([np.broadcast_to(p, (depth, *p.shape[1:])) if p.shape[0] != depth else p for p in planes])
        small = int(stack.max()) < 65536
        tifffile.imwrite(folder / f"{base}_{name}.tif", stack.astype(np.uint16 if small else np.uint32), compression="zlib",
                         imagej=small, resolution=(1 / cal.pixel_size, 1 / cal.pixel_size),
                         metadata={"axes": "TZYX", "spacing": cal.z_spacing, "unit": "micron" if cal.unit == "µm" else cal.unit}
                         if small else None)


# ---------------------------------------------------------------- overlays at full resolution
def full_plane(r: FrameResult, lab: np.ndarray, z: int | None, height: int, width: int) -> np.ndarray:
    """A label image of the analysis as one plane of the image: slice z when the labels are 3D and a slice is
    asked for, otherwise the maximum over Z."""
    small = lab[0] if lab.shape[0] == 1 else (lab[z] if z is not None and z < lab.shape[0] else lab.max(0))
    if r.binning > 1:
        small = np.repeat(np.repeat(small, r.binning, 0), r.binning, 1)
    out = np.zeros((height, width), small.dtype)
    oy, ox = r.offset[1] + r.border, r.offset[0] + r.border
    h, w = min(small.shape[0], height - oy), min(small.shape[1], width - ox)
    out[oy:oy + h, ox:ox + w] = small[:h, :w]
    return out


def edges(lab: np.ndarray) -> np.ndarray:
    """Pixels of a label image that border another label (the outline, one pixel wide, inside the object)."""
    e = np.zeros(lab.shape, bool)
    e[:-1] |= lab[:-1] != lab[1:]
    e[1:] |= lab[1:] != lab[:-1]
    e[:, :-1] |= lab[:, :-1] != lab[:, 1:]
    e[:, 1:] |= lab[:, 1:] != lab[:, :-1]
    return e & (lab > 0)


def render_rgb(planes: list[np.ndarray], display: list[dict], overlays: list[dict]) -> np.ndarray:
    """An RGB picture (uint8) as the viewer shows it: channels added in their colours with their contrast,
    then overlays: {"labels", "colour": (r, g, b), "fill": 0..1, "outline": True, "halo": False}."""
    h, w = planes[0].shape
    rgb = np.zeros((h, w, 3), np.float32)
    for plane, d in zip(planes, display, strict=True):
        if not d.get("on", True):
            continue
        lo, hi = float(d["lo"]), float(d["hi"])
        v = np.clip((plane.astype(np.float32) - lo) / max(hi - lo, 1.0), 0, 1)
        rgb += v[..., None] * (np.asarray(d["colour"], np.float32) / 255)
    rgb = np.clip(rgb, 0, 1)
    for o in overlays:
        lab, colour = o["labels"], np.asarray(o["colour"], np.float32) / 255
        if o.get("fill"):
            m = lab > 0
            rgb[m] = rgb[m] * (1 - o["fill"]) + colour * o["fill"]
        if o.get("outline", True):
            e = edges(lab)
            if o.get("halo"):  # a dark edge around the line, so it shows on any colour
                halo = ndi.binary_dilation(e, iterations=1) & ~e
                rgb[halo] = 0
            rgb[e] = colour
    return (rgb * 255 + 0.5).astype(np.uint8)


# ---------------------------------------------------------------- ImageJ TIFF with overlays
def _contours(lab2d: np.ndarray, ids: np.ndarray | None = None):
    """Per label: the longest outline as (n, 2) x, y points in pixel centres of this plane."""
    for i, sl in enumerate(ndi.find_objects(lab2d), start=1):
        if sl is None or (ids is not None and i not in ids):
            continue
        y0, x0 = sl[0].start - 1, sl[1].start - 1
        mask = np.pad(lab2d[sl] == i, 1)
        cs = measure.find_contours(mask.astype(np.uint8), 0.5)
        if cs:
            c = max(cs, key=len)
            yield i, np.column_stack([c[:, 1] + x0, c[:, 0] + y0])


def frame_rois(r: FrameResult, t: int, frames: int, slices: int, colocalizing: dict | None = None) -> list:
    """ImageJ overlay outlines of one frame: the cells (white) and the foci per channel (their colour).

    2D labels (2D segmentation, or foci found on a projection) are shown in every slice; 3D labels slice by
    slice. Positions are 1-based; 0 means all.
    """
    from roifile import ROI_OPTIONS, ImagejRoi
    no_labels = ~(ROI_OPTIONS.OVERLAY_LABELS | ROI_OPTIONS.OVERLAY_NAMES | ROI_OPTIONS.OVERLAY_BACKGROUNDS
                  | ROI_OPTIONS.SHOW_LABELS)
    out = []
    tpos = t + 1 if frames > 1 else 0

    def add(lab, name, colour):
        for z in range(lab.shape[0]):
            zpos = z + 1 if lab.shape[0] > 1 and slices > 1 else 0
            for i, pts in _contours(lab[z]):
                pts = (pts + 0.5) * r.binning + r.border + np.asarray(r.offset)  # to image pixels (corners)
                roi = ImagejRoi.frompoints(pts, name=name(i))
                roi.position, roi.c_position, roi.z_position, roi.t_position = 0, 0, zpos, tpos  # 0 = all
                roi.stroke_color = bytes((255, *colour))  # no ROI group: ImageJ would colour by group instead
                roi.options &= no_labels  # no number or name drawn at every outline
                out.append(roi)

    add(r.labels, lambda i: f"cell {i}" + (f" t{t + 1}" if frames > 1 else ""), (255, 255, 255))
    for k, ch in enumerate(r.channels):
        if ch.labels is None:
            continue
        with_ = colocalizing.get(ch.channel, set()) if colocalizing else set()
        add(ch.labels, lambda i, c=ch.channel, w=with_: f"ch{c} focus {i}" + (" coloc" if i in w else ""),
            OUTLINE[k % len(OUTLINE)])
    return out


def colocalizing_foci(r: FrameResult) -> dict[int, set]:
    cols = [c for c in r.foci.columns if c.startswith("with_ch")]
    out = {}
    if cols and "channel" in r.foci:
        for ch in r.channels:
            rows = r.foci[r.foci["channel"] == ch.channel]
            hit = rows[cols].fillna(0).astype(float).to_numpy().any(axis=1)
            out[ch.channel] = set(rows["focus"].to_numpy()[hit].tolist())
    return out


def write_imagej(path: Path, data: np.ndarray, info, frames: dict[int, FrameResult], colours: list, ranges=None) -> int:
    """The image (T, C, Z, Y, X) as an ImageJ hyperstack TIFF with the outlines of the analysis as overlay
    (only for the frames in `frames`) and the channel colours. Returns the number of outlines.

    Not compressed: ImageJ and Fiji drop the overlay of a compressed TIFF (tested with Fiji, 2026-11-17).
    """
    T, Z = data.shape[0], data.shape[2]
    rois = []
    for t, r in sorted(frames.items()):
        rois += frame_rois(r, t, T, Z, colocalizing_foci(r))
    stack = np.ascontiguousarray(data.transpose(0, 2, 1, 3, 4))  # ImageJ order: T, Z, C, Y, X
    if stack.dtype not in (np.uint8, np.uint16, np.float32):
        stack = stack.astype(np.float32)
    luts = []
    for col in colours:
        ramp = np.arange(256, dtype=np.float64) / 255
        luts.append(np.stack([(ramp * c).astype(np.uint8) for c in col]))
    meta = {"axes": "TZCYX", "mode": "composite", "LUTs": luts, "unit": "micron" if info.unit == "µm" else info.unit}
    if Z > 1 and info.z_spacing:
        meta["spacing"] = info.z_spacing
    if T > 1 and info.frame_interval:
        meta["finterval"] = info.frame_interval
    if ranges:
        meta["Ranges"] = tuple(float(v) for lohi in ranges for v in lohi)
    if rois:
        meta["Overlays"] = [roi.tobytes() for roi in rois]
    path.parent.mkdir(parents=True, exist_ok=True)
    px = info.pixel_size or 1.0
    tmp = path.with_suffix(".tmp.tif")
    tifffile.imwrite(tmp, stack, imagej=True, resolution=(1 / px, 1 / px), metadata=meta)
    tmp.replace(path)
    return len(rois)
