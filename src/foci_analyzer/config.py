"""Settings: one definition for the form in the window, validation and the TOML files.

Every setting is a dataclass field made with `opt()`. Its metadata (label, help, advanced, choices, limits,
unit) is all the window needs to draw the form, so a new setting is added in one place only.
"""

from __future__ import annotations

import tomllib
from dataclasses import MISSING, dataclass, field, fields
from pathlib import Path
from typing import Any

import tomli_w


def opt(default: Any, label: str, help: str = "", *, advanced: bool = False, choices: list[str] | None = None,
        min: float | None = None, max: float | None = None, step: float | None = None, unit: str = "",
        slider: bool = False):
    """A setting with what the window needs to show it. `help` is the context help (the ? in the form)."""
    meta = {"label": label, "help": help, "advanced": advanced, "choices": choices, "min": min, "max": max,
            "step": step, "unit": unit, "slider": slider}
    return field(default=default, metadata=meta)


@dataclass
class ImageSettings:
    """Image"""
    nuclei_channel: int = opt(1, "Nuclei channel", "The channel with the nuclear stain (for example DAPI). "
                              "Foci are counted per nucleus. 0 = no nuclei channel.", min=0)
    cell_channel: int = opt(0, "Cell channel", "A channel that shows whole cells (cytoplasm or membrane). When "
                            "set, cells are segmented instead of nuclei and foci are counted per cell. "
                            "0 = none. With both channels at 0 the whole image is treated as one cell.", min=0)
    xy_binning: int = opt(1, "XY binning", "Average n × n pixels before the analysis. Use it for images "
                          "sampled much finer than the foci size; it reduces noise and time.", advanced=True,
                          min=1)
    crop_border: int = opt(0, "Remove border", "Pixels removed from every image edge before the analysis, for "
                           "example a dark stitching border.", advanced=True, min=0, unit="px")


@dataclass
class SegmentationSettings:
    """Nuclei / cells"""
    method: str = opt("threshold", "Method", "How nuclei or cells are found. Cellpose, StarDist and InstanSeg "
                      "are neural networks and need the models to be installed. 'threshold' smooths the image, "
                      "applies an Otsu threshold and splits touching objects on the distance map. 'labels from "
                      "file' reads label images made elsewhere.",
                      choices=["threshold", "cellpose", "stardist", "instanseg", "labels from file"])
    dimension: str = opt("2D", "Segment in", "2D segments the maximum projection of the stack; the outlines "
                         "then hold for all slices. This is faster and often more robust, and can be combined "
                         "with foci detection in 3D. 3D segments the stack itself.", choices=["2D", "3D"])
    min_diameter: float = opt(4.0, "Smallest diameter", "Objects with a smaller diameter (of a circle with "
                              "the same area) are removed.", min=0, unit="µm")
    max_diameter: float = opt(50.0, "Largest diameter", "Objects with a larger diameter are removed. "
                              "0 = no limit.", min=0, unit="µm")
    exclude_on_edges: bool = opt(True, "Remove objects on the image edge", "A nucleus that is cut by the image "
                                 "edge has an incomplete foci count.", advanced=True)
    model: str = opt("", "Model", "Model name or path; empty = the default of the method.", advanced=True)
    diameter: float = opt(0.0, "Expected diameter", "For Cellpose: the typical diameter. 0 = automatic.",
                          advanced=True, min=0, unit="µm")
    probability: float = opt(0.5, "Probability threshold", "Lower finds more and larger objects.",
                             advanced=True, min=0, max=1, step=0.05)
    labels_folder: str = opt("", "Labels folder", "For 'labels from file': a folder with label images (TIFF) "
                             "that have the same name as the images.", advanced=True)
    python: str = opt("", "Python of the models", "Path of the python.exe of another environment that holds the "
                      "model (for example an existing Cellpose environment). Empty = the environment of "
                      "Foci Analyzer.", advanced=True)


@dataclass
class FociGeneralSettings:
    """Foci"""
    dimension: str = opt("3D", "Detect foci in", "For stacks. 3D detects foci in the stack with filters that "
                         "work in X, Y and Z. A projection detects them in one plane, which is often more "
                         "reliable for thin samples and few slices. Ignored for 2D images.",
                         choices=["3D", "max projection", "extended depth of focus", "sum projection"])
    outside_distance: float = opt(0.0, "Distance outside cells", "Foci up to this distance outside a nucleus "
                                  "or cell still count for it. -1 = every focus counts for the nearest cell.",
                                  advanced=True, min=-1, unit="µm")
    z_elongation: float = opt(2.0, "Axial elongation", "How much longer a focus is in Z than in XY because "
                              "of the microscope optics. Used for the filter size in Z (3D detection).",
                              advanced=True, min=0.5, max=10, step=0.5, unit="×")


@dataclass
class FociChannelSettings:
    """Foci channel"""
    channel: int = opt(2, "Channel", "The channel with the foci.", min=1)
    method: str = opt("watershed", "Method", "watershed: the image is filtered so that spots of the foci size "
                      "stand out, local maxima above the threshold become foci, and a seeded watershed gives "
                      "every focus its extent. spotiflow: a neural network finds the spots and gives each a "
                      "probability; the extents are found in the same way.",
                      choices=["watershed", "spotiflow"])
    size: float = opt(0.5, "Foci size", "Expected diameter of a focus (full width at half maximum). It sets "
                      "the filter scale and the smallest distance between two foci.", min=0.05, step=0.05,
                      unit="µm")
    sensitivity: float = opt(8.0, "Threshold", "A local maximum is a focus when it is this many times higher "
                             "than the noise of the filtered image. The noise is 1.4826 × the median absolute "
                             "deviation of the filtered values inside the cells, which is hardly influenced "
                             "by the foci themselves. Lower finds more foci. The viewer follows this slider "
                             "directly.", min=1, max=50, step=0.1, unit="× noise", slider=True)
    probability: float = opt(0.5, "Spotiflow probability", "Spots below this probability are left out.",
                             advanced=True, min=0, max=1, step=0.01, slider=True)
    filter: str = opt("DoG", "Filter", "DoG: difference of two Gaussians around the foci size. LoG: Laplacian "
                      "of Gaussian at the foci size; slightly sharper, more sensitive to noise. wavelet: "
                      "product of two à trous wavelet planes; good for foci of mixed size on uneven "
                      "background. All work in 3D for 3D detection.", advanced=True,
                      choices=["DoG", "LoG", "wavelet"])
    extent: str = opt("half maximum", "Focus extent", "Where a focus ends: at half of its own peak height "
                      "(size does not depend on the threshold), or at the threshold.", advanced=True,
                      choices=["half maximum", "threshold"])
    min_size: int = opt(2, "Smallest focus", "Smaller foci are removed.", advanced=True, min=1, unit="voxels")
    max_size: int = opt(9999, "Largest focus", "Larger foci are removed.", advanced=True, min=1, unit="voxels")
    threshold_scope: str = opt("per image", "Noise from", "per image: one noise value from all cells of the "
                               "image. per cell: every cell its own noise. fixed value: the same absolute "
                               "threshold for all images, which is the fair choice when conditions differ "
                               "strongly in foci numbers.", advanced=True,
                               choices=["per image", "per cell", "fixed value"])
    fixed_threshold: float = opt(0.0, "Fixed threshold", "For 'fixed value': the threshold in units of the "
                                 "filtered image. The viewer shows the value that the current image uses.",
                                 advanced=True, min=0)
    spotiflow_model: str = opt("general", "Spotiflow model", "Name of a pretrained model or a path.",
                               advanced=True)


@dataclass
class ColocalizationSettings:
    """Colocalization"""
    enabled: bool = opt(True, "Colocalize foci channels", "Used when there are two or more foci channels; "
                        "every pair of channels is compared.")
    method: str = opt("overlap", "Method", "overlap: two foci colocalize when they share voxels. centroid "
                      "distance: when their centres are closer than the largest distance.",
                      choices=["overlap", "centroid distance"])
    min_overlap: int = opt(1, "Smallest overlap", advanced=True, min=1, unit="voxels")
    max_distance: float = opt(0.3, "Largest distance", advanced=True, min=0, step=0.05, unit="µm")
    chance_repeats: int = opt(0, "Chance control repeats", "The colocalization expected by chance: the foci of "
                              "the second channel are placed at random positions in the cell this many times "
                              "(0 = off).", advanced=True, min=0)


@dataclass
class TimelapseSettings:
    """Time-lapse"""
    tracking: str = opt("off", "Track cells", "Follow nuclei or cells over the frames. overlap: a cell is "
                        "linked to the cell it overlaps most in the previous frame. trackastra: a neural "
                        "network that also handles divisions.", choices=["off", "overlap", "trackastra"])


@dataclass
class Settings:
    image: ImageSettings = field(default_factory=ImageSettings)
    segmentation: SegmentationSettings = field(default_factory=SegmentationSettings)
    foci_general: FociGeneralSettings = field(default_factory=FociGeneralSettings)
    foci: list[FociChannelSettings] = field(default_factory=lambda: [FociChannelSettings()])
    colocalization: ColocalizationSettings = field(default_factory=ColocalizationSettings)
    timelapse: TimelapseSettings = field(default_factory=TimelapseSettings)


REPEATED = {"foci": FociChannelSettings}  # sections that are a list of the same group


def _field_schema(f) -> dict:
    default = f.default if f.default is not MISSING else None
    kind = {bool: "bool", int: "int", float: "float", str: "str"}[type(default)]
    if f.metadata.get("choices"):
        kind = "choice"
    return {"key": f.name, "type": kind, "default": default, **f.metadata}


def schema() -> list[dict]:
    """The sections and their settings, in the order of the form."""
    out = []
    for f in fields(Settings):
        cls = REPEATED.get(f.name) or f.default_factory
        out.append({"key": f.name, "title": (cls.__doc__ or f.name).strip(), "repeat": f.name in REPEATED,
                    "fields": [_field_schema(g) for g in fields(cls)]})
    return out


def _coerce(f, value):
    """The value as the type of the setting, inside its limits; ValueError when it cannot be used."""
    default, m = f.default, f.metadata
    if isinstance(default, bool):
        if not isinstance(value, bool):
            raise ValueError(f"{f.name}: expected true or false")
        return value
    if isinstance(default, (int, float)):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{f.name}: expected a number")
        value = type(default)(value)
        if m["min"] is not None and value < m["min"]:
            raise ValueError(f"{f.name}: smaller than {m['min']}")
        if m["max"] is not None and value > m["max"]:
            raise ValueError(f"{f.name}: larger than {m['max']}")
        return value
    if not isinstance(value, str):
        raise ValueError(f"{f.name}: expected text")
    if m["choices"] and value not in m["choices"]:
        raise ValueError(f"{f.name}: '{value}' is not one of {m['choices']}")
    return value


def _group_from_dict(cls, data: dict):
    known = {f.name: f for f in fields(cls)}
    return cls(**{k: _coerce(known[k], v) for k, v in (data or {}).items() if k in known})


def from_dict(data: dict) -> Settings:
    """Settings from a dictionary; missing settings get their default, unknown ones are ignored."""
    s = Settings()
    for f in fields(Settings):
        if f.name not in data:
            continue
        if f.name in REPEATED:
            groups = [_group_from_dict(REPEATED[f.name], d) for d in data[f.name]]
            if not groups:
                raise ValueError(f"{f.name}: at least one is needed")
            setattr(s, f.name, groups)
        else:
            setattr(s, f.name, _group_from_dict(f.default_factory, data[f.name]))
    return s


def to_dict(s: Settings) -> dict:
    from dataclasses import asdict
    return asdict(s)


def load(path: Path) -> Settings:
    return from_dict(tomllib.loads(Path(path).read_text(encoding="utf-8")))


def save(s: Settings, path: Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(tomli_w.dumps(to_dict(s)), encoding="utf-8")
