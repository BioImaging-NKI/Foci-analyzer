"""Reading images. Every image becomes an array in the order T, C, Z, Y, X with its calibration.

TIFF (also ImageJ and OME), Zeiss CZI (every scene is an image) and Leica LIF (every series is an image).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import tifffile


def _czi():
    try:
        from pylibCZIrw import czi
    except ImportError as e:
        raise ValueError("CZI files cannot be read on this computer (the CZI reader pylibCZIrw has no version for "
                         "it, for example an Intel Mac); export the images as TIFF") from e
    return czi

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
    colors: list | None = None  # per channel the color the file gives it ("#rrggbb"), None for gray or unknown

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


# ---- channel colors as the microscope or Fiji stored them; gray counts as no color
_NAMED = {"red": (255, 0, 0), "green": (0, 255, 0), "blue": (0, 0, 255), "cyan": (0, 255, 255),
          "magenta": (255, 0, 255), "yellow": (255, 255, 0)}


def _hex(rgb) -> str | None:
    r, g, b = (int(v) for v in rgb)
    return None if r == g == b else f"#{r:02x}{g:02x}{b:02x}"


def _colors(found: list, channels: int) -> list | None:
    out = [found[c] if c < len(found) else None for c in range(channels)]
    return out if any(out) else None


def _tiff_colors(tf, channels: int) -> list | None:
    try:
        luts = (tf.imagej_metadata or {}).get("LUTs")
        if luts is not None:
            return _colors([_hex(np.asarray(lut)[:, -1]) for lut in luts], channels)
        if tf.is_ome:
            import re
            found = []
            for m in re.finditer(r"<(?:\w+:)?Channel\b[^>]*>", tf.ome_metadata or ""):
                c = re.search(r'\bColor="(-?\d+)"', m.group(0))
                v = int(c.group(1)) & 0xFFFFFFFF if c else None  # RGBA as a signed 32-bit number
                found.append(_hex(((v >> 24) & 255, (v >> 16) & 255, (v >> 8) & 255)) if v is not None else None)
            return _colors(found, channels)
    except (ValueError, TypeError, KeyError, IndexError, AttributeError):  # colors never stop an image opening
        return None
    return None


def list_series(path: str | Path) -> list[str]:
    """The names of the images in a file; one empty name for a file with a single image."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".lif":
        from readlif.reader import LifFile
        return [im["name"] for im in LifFile(str(path)).image_list]
    if ext == ".czi":
        czi = _czi()
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
        info = _info(path, path.stem, arr, px, ij.get("spacing"), 0, ij.get("finterval"))
        info.colors = _tiff_colors(tf, info.channels)
        return arr, info


def _read_czi(path: Path, series: int):
    czi = _czi()
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
        try:  # "#AARRGGBB" per channel
            chans = doc.metadata["ImageDocument"]["Metadata"]["DisplaySetting"]["Channels"]["Channel"]
            chans = chans if isinstance(chans, list) else [chans]
            found = [_hex(int(str(c.get("Color"))[-6:][i:i + 2], 16) for i in (0, 2, 4)) if c.get("Color") else None
                     for c in chans]
        except (KeyError, TypeError, ValueError, AttributeError):
            found = []
    name = path.stem + (f" scene {series + 1}" if len(scenes) > 1 else "")
    info = _info(path, name, arr, scale.get("X"), scale.get("Z"), series)
    info.colors = _colors(found, info.channels)
    return arr, info


def _read_lif(path: Path, series: int):
    from readlif.reader import LifFile
    lif = LifFile(str(path))
    im = lif.get_image(series)
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
    info = _info(path, f"{path.stem} {im.name}", arr, 1 / sx if sx else None, 1 / abs(sz) if sz else None,
                 series, 1 / st if st else None)
    try:  # the LUT names of the channels of this image (the n-th image of that name in the file)
        same = [e for e in lif.xml_root.iter("Element") if e.find("Data/Image") is not None and e.get("Name") == im.name]
        nth = [x["name"] for x in lif.image_list[:series]].count(im.name)
        found = [_hex(_NAMED[n]) if (n := ch.get("LUTName", "").lower()) in _NAMED else None
                 for ch in same[nth].iter("ChannelDescription")]
        info.colors = _colors(found, info.channels)
    except (IndexError, AttributeError):
        pass
    return arr, info


def read_image(path: str | Path, series: int = 0) -> tuple[np.ndarray, ImageInfo]:
    """The pixels (T, C, Z, Y, X) and the description of an image file."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".czi":
        return _read_czi(path, series)
    if ext == ".lif":
        return _read_lif(path, series)
    return _read_tiff(path)
