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
    if channels > 1 and len(set(out)) == 1:  # every channel the same color (some LIF files: all blue): not a choice
        return None
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


_LENGTH = {"µm": 1.0, "um": 1.0, "nm": 1e-3, "mm": 1e3, "cm": 1e4, "m": 1e6, "Å": 1e-4}
_TIME = {"s": 1.0, "ms": 1e-3, "µs": 1e-6, "us": 1e-6, "min": 60.0, "h": 3600.0}


def _ome_calibration(xml: str) -> tuple[float | None, float | None, float | None]:
    """Pixel size and Z spacing in µm and the frame interval in seconds, from the first image of OME-XML."""
    import re
    m = re.search(r"<(?:\w+:)?Pixels\b[^>]*>", xml)
    if not m:
        return None, None, None
    tag = m.group(0)

    def value(name, units, default):
        v = re.search(rf'\b{name}="([^"]+)"', tag)
        u = re.search(rf'\b{name}Unit="([^"]+)"', tag)
        try:
            return float(v.group(1)) * units[u.group(1) if u else default] if v else None
        except (ValueError, KeyError):
            return None
    return value("PhysicalSizeX", _LENGTH, "µm"), value("PhysicalSizeZ", _LENGTH, "µm"), value("TimeIncrement", _TIME, "s")


def list_series(path: str | Path) -> list[str]:
    """The names of the images in a file; one empty name for a file with a single image."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".lif":
        import liffile
        with liffile.LifFile(path) as lif:
            return [im.path for im in _lif_images(lif)]
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
        z, interval = ij.get("spacing"), ij.get("finterval")
        if tf.is_ome and not ij:  # an OME-TIFF without ImageJ's description: the calibration from the OME data
            px, z, interval = _ome_calibration(tf.ome_metadata or "")
        info = _info(path, path.stem, arr, px, z, 0, interval)
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
        k = first.shape[-1] if first.ndim == 3 else 1  # 3 for a color (RGB) image: each color becomes a channel
        arr = np.zeros((nt, nc * k, nz, rect.h, rect.w), first.dtype)
        for t in range(nt):
            for c in range(nc):
                for z in range(nz):
                    plane = doc.read(plane={"T": box["T"][0] + t, "C": box["C"][0] + c, "Z": box["Z"][0] + z},
                                     scene=scene, roi=(rect.x, rect.y, rect.w, rect.h))
                    plane = plane.reshape(rect.h, rect.w, -1)
                    for j in range(k):  # the library gives the colors as blue, green, red
                        arr[t, c * k + j, z] = plane[..., k - 1 - j]
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
        if k == 3:  # a color image: its own red, green and blue
            found = ["#ff0000", "#00ff00", "#0000ff"] * nc
    name = path.stem + (f" scene {series + 1}" if len(scenes) > 1 else "")
    info = _info(path, name, arr, scale.get("X"), scale.get("Z"), series)
    info.colors = _colors(found, info.channels)
    return arr, info


def tile_offsets(tiles: list[dict], step: float, w: int, h: int, flags: dict | None = None) -> np.ndarray:
    """(n, 2) pixel offset (x, y) of every tile of a Leica tile scan, from its stage position.

    tiles: the attributes of the Leica <Tile> elements (FieldX, FieldY: place in the grid; PosX, PosY: stage position
    in m). step: the size of a pixel in µm. flags: Leica's FlipX, FlipY and SwapXY (True or False) of the scan: X and Y
    swapped, then X and/or Y negated. The image then runs along X and Y (checked against Fiji's stitching of a Leica
    tile scan with all flags off). Without stage positions the tiles go on the grid, without overlap."""
    f = flags or {}
    try:
        pos = np.array([[float(t["PosX"]), float(t["PosY"])] for t in tiles]) * 1e6 / step
    except (KeyError, TypeError, ValueError):
        pos = None
    if pos is None or not np.isfinite(pos).all():
        out = np.array([[float(t.get("FieldX", 0)) * w, float(t.get("FieldY", 0)) * h] for t in tiles])
    else:
        if f.get("SwapXY"):
            pos = pos[:, ::-1]
        out = pos * (-1 if f.get("FlipX") else 1, -1 if f.get("FlipY") else 1)
    return np.rint(out - out.min(axis=0)).astype(int)


def _lif_images(lif) -> list:
    """The images of a .lif file, in file order, without the raw data that are not plain images (FLIM decays)."""
    return [im for im in lif.images if set(im.dims) <= set("TMCZYX")]


def _read_lif(path: Path, series: int):
    import liffile
    with liffile.LifFile(path) as lif:
        im = _lif_images(lif)[series]
        data = np.asarray(im.asarray())
        coords, scan, name, xml, axes = im.coords, im.tilescan, im.path, im.xml_element, "".join(im.dims)
    for a in "TMCZ":
        if a not in axes:
            data, axes = data[np.newaxis], a + axes
    data = np.transpose(data, [axes.index(a) for a in "TMCZYX"])
    nt, nm, nc, nz, h, w = data.shape
    step = {a: abs(float(np.mean(np.diff(v)))) for a, v in coords.items() if a in "XZT" and len(v) > 1}
    px = step["X"] * 1e6 if step.get("X") else None
    if nm > 1:  # a tile scan: every tile where the stage was (a later tile covers the overlap of an earlier one)
        if scan is not None and len(scan.tiles) == nm:
            tiles = [{"FieldX": t["field_x"], "FieldY": t["field_y"], "PosX": t["pos_x"], "PosY": t["pos_y"]}
                     for t in scan.tiles]
            flags = {"FlipX": scan.flip_x, "FlipY": scan.flip_y, "SwapXY": scan.swap_xy}
        else:
            tiles, flags = [{"FieldX": i, "FieldY": 0} for i in range(nm)], None
        offsets = tile_offsets(tiles, px or 1.0, w, h, flags)
        arr = np.zeros((nt, nc, nz, int(offsets[:, 1].max()) + h, int(offsets[:, 0].max()) + w), data.dtype)
        for m, (ox, oy) in enumerate(offsets):
            arr[..., oy:oy + h, ox:ox + w] = data[:, m]
    else:
        arr = data[:, 0]
    info = _info(path, f"{path.stem} {name}", arr, px, step["Z"] * 1e6 if step.get("Z") else None, series,
                 step.get("T"))
    try:  # the LUT names of the channels of this image
        found = [_hex(_NAMED[lut]) if (lut := ch.get("LUTName", "").lower()) in _NAMED else None
                 for ch in xml.iter("ChannelDescription")]
        info.colors = _colors(found, info.channels)
    except AttributeError:
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
