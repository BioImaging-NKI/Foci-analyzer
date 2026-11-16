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
        min: float | None = None, max: float | None = None, step: float | None = None, unit: str = ""):
    """A setting with what the window needs to show it."""
    meta = {"label": label, "help": help, "advanced": advanced, "choices": choices, "min": min, "max": max,
            "step": step, "unit": unit}
    return field(default=default, metadata=meta)


@dataclass
class ImageSettings:
    """Image"""
    nuclei_channel: int = opt(1, "Nuclei channel", "Channel with the nuclei (0 = none).", min=0)
    cell_channel: int = opt(0, "Cell channel", "Channel with cytoplasm or membrane (0 = none).", min=0)
    xy_binning: int = opt(1, "XY binning", "Average n × n pixels before the analysis.", advanced=True, min=1)
    crop_border: int = opt(0, "Remove border", "Pixels removed from every image edge.", advanced=True, min=0,
                           unit="px")


@dataclass
class SegmentationSettings:
    """Nuclei / cells"""
    method: str = opt("cellpose", "Method", "How nuclei or cells are found.",
                      choices=["cellpose", "stardist", "instanseg", "threshold", "labels from file"])
    dimension: str = opt("2D", "Segment in", "2D uses a projection of the stack; the labels then hold for "
                         "all slices.", choices=["2D", "3D"])
    min_diameter: float = opt(4.0, "Smallest diameter", "Smaller objects are removed.", min=0, unit="µm")
    max_diameter: float = opt(50.0, "Largest diameter", "Larger objects are removed.", min=0, unit="µm")
    exclude_on_edges: bool = opt(True, "Remove objects on the image edge", advanced=True)
    model: str = opt("", "Model", "Model name or path; empty = the default of the method.", advanced=True)
    diameter: float = opt(0.0, "Expected diameter", "0 = automatic.", advanced=True, min=0, unit="µm")
    probability: float = opt(0.5, "Probability / flow threshold", advanced=True, min=0, max=1, step=0.05)
    labels_folder: str = opt("", "Labels folder", "For 'labels from file': label images with the same name "
                             "as the images.", advanced=True)


@dataclass
class FociGeneralSettings:
    """Foci"""
    dimension: str = opt("3D", "Detect foci in", "For stacks: in 3D, on a projection, or in every slice "
                         "separately. Ignored for 2D images.",
                         choices=["3D", "max projection", "extended depth of focus", "sum projection",
                                  "every slice"])
    outside_distance: float = opt(0.0, "Distance outside cells", "Foci up to this distance outside a cell "
                                  "still count (-1 = whole image).", advanced=True, min=-1, unit="µm")


@dataclass
class FociChannelSettings:
    """Foci channel"""
    channel: int = opt(2, "Channel", min=1)
    method: str = opt("watershed", "Method", "Seeded watershed, or the Spotiflow network.",
                      choices=["watershed", "spotiflow"])
    size: float = opt(0.5, "Foci size", "Expected diameter of a focus.", min=0.05, step=0.05, unit="µm")
    sensitivity: float = opt(5.0, "Threshold", "Watershed: times the noise of the filtered image. Lower finds "
                             "more foci.", min=1, max=50, step=0.1, unit="× noise")
    probability: float = opt(0.5, "Spotiflow probability", "Spots below this probability are left out.",
                             advanced=True, min=0, max=1, step=0.01)
    min_size: int = opt(3, "Smallest focus", advanced=True, min=1, unit="voxels")
    max_size: int = opt(9999, "Largest focus", advanced=True, min=1, unit="voxels")
    extent: str = opt("half maximum", "Focus extent", "Where a focus ends: at half of its own peak, or at the "
                      "threshold.", advanced=True, choices=["half maximum", "threshold"])
    threshold_scope: str = opt("per image", "Threshold from", "A fixed value is the fair choice when "
                               "conditions are compared.", advanced=True,
                               choices=["per image", "per cell", "fixed value"])
    fixed_threshold: float = opt(0.0, "Fixed threshold", "In units of the filtered image.", advanced=True, min=0)
    spotiflow_model: str = opt("general", "Spotiflow model", "Name or path.", advanced=True)


@dataclass
class ColocalizationSettings:
    """Colocalization"""
    enabled: bool = opt(True, "Colocalize foci channels", "Used when there are two or more foci channels.")
    method: str = opt("overlap", "Method", choices=["overlap", "centroid distance"])
    min_overlap: int = opt(1, "Smallest overlap", advanced=True, min=1, unit="voxels")
    max_distance: float = opt(0.3, "Largest distance", advanced=True, min=0, step=0.05, unit="µm")
    chance_repeats: int = opt(0, "Chance control repeats", "Random placements per cell (0 = off).",
                              advanced=True, min=0)


@dataclass
class TimelapseSettings:
    """Time-lapse"""
    tracking: str = opt("off", "Track cells", choices=["off", "overlap", "trackastra"])


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
