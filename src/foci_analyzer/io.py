"""Reading images. Every image becomes an array in the order T, C, Z, Y, X with its calibration."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import tifffile

EXTENSIONS = (".tif", ".tiff")  # more formats come through bioio (plan.md, section 3)
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


def read_image(path: str | Path) -> tuple[np.ndarray, ImageInfo]:
    """The pixels (T, C, Z, Y, X) and the description of a TIFF file."""
    path = Path(path)
    with tifffile.TiffFile(path) as tf:
        series = tf.series[0]
        arr = _to_tczyx(series.asarray(), series.axes)
        ij = tf.imagej_metadata or {}
        tags = tf.pages[0].tags
        pixel_size = None
        if "XResolution" in tags:
            num, den = tags["XResolution"].value
            if num and den and num != den:  # 1/1 is what programs write when there is no calibration
                pixel_size = den / num
        unit = str(ij.get("unit", "")).strip().lower()
        calibrated = pixel_size is not None and unit in _MICRON
        z = ij.get("spacing")
    t, c, nz, h, w = arr.shape
    info = ImageInfo(path=str(path), name=path.stem, frames=t, channels=c, slices=nz, height=h, width=w,
                     dtype=str(arr.dtype), pixel_size=pixel_size if calibrated else None,
                     z_spacing=float(z) if (calibrated and z and nz > 1) else None,
                     unit="µm" if calibrated else "px")
    return arr, info
