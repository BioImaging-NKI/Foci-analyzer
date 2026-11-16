"""Reading images. Every image becomes an array in the order T, C, Z, Y, X with its calibration.

TIFF (also ImageJ and OME), Zeiss CZI (every scene is an image) and Leica LIF (every series is an image).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import tifffile

EXTENSIONS = (".tif", ".tiff", ".czi", ".lif")
_MICRON = {"micron", "microns", "um", "µm", "μm", "\\u00b5m"}
_AXES = "TCZYX"


@dataclass
class ImageInfo:
    path: str
    name: str
    frames: int
    channels: int
    slices: int
    height: int
    width: int
    dtype: str
    pixel_size: float | None  # XY, in `unit`; None when the file has no calibration
    z_spacing: float | None
    unit: str  # "µm" or "px"
    series: int = 0
    frame_interval: float | None = None  # seconds

    def to_dict(self) -> dict:
        return asdict(self)


def _to_tczyx(arr: np.ndarray, axes: str) -> np.ndarray:
    # S (samples of an RGB image) counts as channels; unnamed stack axes (I, Q) as slices
    axes = axes.replace("S", "C").replace("I", "Z").replace("Q", "Z")
    if len(set(axes)) != len(axes) or not set(axes) <= set(_AXES):
        raise ValueError(f"cannot interpret the image axes '{axes}'")
    for a in _AXES:
        if a not in axes:
            arr, axes = arr[np.newaxis], a + axes
    return np.transpose(arr, [axes.index(a) for a in _AXES])


def _info(path: Path, name: str, arr: np.ndarray, px, z, series=0, interval=None) -> ImageInfo:
    t, c, nz, h, w = arr.shape
    return ImageInfo(path=str(path), name=name, frames=t, channels=c, slices=nz, height=h, width=w,
                     dtype=str(arr.dtype), pixel_size=px or None, z_spacing=float(z) if (px and z and nz > 1) else None,
                     unit="µm" if px else "px", series=series, frame_interval=interval)


def list_series(path: str | Path) -> list[str]:
    """The names of the images in a file; one empty name for a file with a single image."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".lif":
        from readlif.reader import LifFile
        return [im["name"] for im in LifFile(str(path)).image_list]
    if ext == ".czi":
        from pylibCZIrw import czi
        with czi.open_czi(str(path)) as doc:
            scenes = sorted(doc.scenes_bounding_rectangle)
        return [f"scene {s + 1}" for s in scenes] if len(scenes) > 1 else [""]
    return [""]


def _read_tiff(path: Path):
    with tifffile.TiffFile(path) as tf:
        series = tf.series[0]
        arr = _to_tczyx(series.asarray(), series.axes)
        ij = tf.imagej_metadata or {}
        tags = tf.pages[0].tags
        px = None
        if "XResolution" in tags:
            num, den = tags["XResolution"].value
            if num and den and num != den:  # 1/1 is what programs write when there is no calibration
                px = den / num
        if str(ij.get("unit", "")).strip().lower() not in _MICRON:
            px = None
        return arr, _info(path, path.stem, arr, px, ij.get("spacing"), 0, ij.get("finterval"))


def _read_czi(path: Path, series: int):
    from pylibCZIrw import czi
    with czi.open_czi(str(path)) as doc:
        box = doc.total_bounding_box
        scenes = sorted(doc.scenes_bounding_rectangle)
        scene = scenes[series] if scenes else None
        rect = doc.scenes_bounding_rectangle[scene] if scenes else doc.total_bounding_rectangle
        nt, nc, nz = (box[a][1] - box[a][0] for a in "TCZ")
        first = doc.read(plane={"T": box["T"][0], "C": box["C"][0], "Z": box["Z"][0]}, scene=scene,
                         roi=(rect.x, rect.y, rect.w, rect.h))
        arr = np.zeros((nt, nc, nz, rect.h, rect.w), first.dtype)
        for t in range(nt):
            for c in range(nc):
                for z in range(nz):
                    plane = doc.read(plane={"T": box["T"][0] + t, "C": box["C"][0] + c, "Z": box["Z"][0] + z},
                                     scene=scene, roi=(rect.x, rect.y, rect.w, rect.h))
                    arr[t, c, z] = plane[..., 0]
        try:
            items = doc.metadata["ImageDocument"]["Metadata"]["Scaling"]["Items"]["Distance"]
            scale = {d["@Id"]: float(d["Value"]) * 1e6 for d in (items if isinstance(items, list) else [items])}
        except (KeyError, TypeError, ValueError):
            scale = {}
    name = path.stem + (f" scene {series + 1}" if len(scenes) > 1 else "")
    return arr, _info(path, name, arr, scale.get("X"), scale.get("Z"), series)


def _read_lif(path: Path, series: int):
    from readlif.reader import LifFile
    im = LifFile(str(path)).get_image(series)
    d = im.dims
    arr = None
    for t in range(d.t):
        for c in range(im.channels):
            for z in range(d.z):
                plane = np.asarray(im.get_frame(z=z, t=t, c=c))
                if arr is None:
                    arr = np.zeros((d.t, im.channels, d.z, *plane.shape), plane.dtype)
                arr[t, c, z] = plane
    sx, _, sz, st = (list(im.scale) + [None] * 4)[:4]
    return arr, _info(path, f"{path.stem} {im.name}", arr, 1 / sx if sx else None, 1 / abs(sz) if sz else None,
                      series, 1 / st if st else None)


def read_image(path: str | Path, series: int = 0) -> tuple[np.ndarray, ImageInfo]:
    """The pixels (T, C, Z, Y, X) and the description of an image file."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".czi":
        return _read_czi(path, series)
    if ext == ".lif":
        return _read_lif(path, series)
    return _read_tiff(path)
