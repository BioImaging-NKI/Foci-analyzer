"""An experiment: conditions, each with one or more folders (or single files) of images."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from fnmatch import fnmatch
from pathlib import Path

from .io import EXTENSIONS, list_series


@dataclass
class Condition:
    name: str
    folders: list[str] = field(default_factory=list)  # folders, or image files
    pattern: str = "*"  # file-name pattern, without the extension
    recursive: bool = False


@dataclass
class ImageEntry:
    path: str
    name: str
    condition: str
    folder: str
    series: int = 0
    id: str = ""  # what the window uses to name the image: the path, with the series for files with several

    def to_dict(self) -> dict:
        return asdict(self)


def conditions_from_dicts(items: list[dict]) -> list[Condition]:
    out, seen = [], set()
    for d in items:
        name = str(d.get("name", "")).strip()
        if not name:
            raise ValueError("a condition needs a name")
        if name in seen:
            raise ValueError(f"two conditions are called '{name}'")
        seen.add(name)
        folders = [str(f).strip().strip('"') for f in d.get("folders", []) if str(f).strip()]
        out.append(Condition(name, folders, str(d.get("pattern") or "*"), bool(d.get("recursive", False))))
    return out


def scan(conditions: list[Condition]) -> tuple[list[ImageEntry], list[str]]:
    """The images of all conditions, and remarks about folders that could not be used."""
    images, remarks, seen = [], [], {}
    for cond in conditions:
        for folder in cond.folders:
            root = Path(folder)
            if root.is_file():
                found = [root] if root.suffix.lower() in EXTENSIONS else []
                if not found:
                    remarks.append(f"{cond.name}: not an image that can be read: {folder}")
            elif not root.is_dir():
                remarks.append(f"{cond.name}: folder not found: {folder}")
                continue
            else:
                files = root.rglob("*") if cond.recursive else root.iterdir()
                found = sorted(p for p in files if p.is_file() and p.suffix.lower() in EXTENSIONS
                               and fnmatch(p.stem, cond.pattern))
                if not found:
                    remarks.append(f"{cond.name}: no images in {folder}")
            for p in found:
                key = str(p.resolve()).lower()
                if key in seen:
                    remarks.append(f"{p.name} is in both '{seen[key]}' and '{cond.name}'; kept in the first")
                    continue
                seen[key] = cond.name
                try:
                    series = list_series(p)
                except Exception as e:  # noqa: BLE001  a file that cannot be opened is a remark, not a failure
                    remarks.append(f"{p.name} cannot be read: {e}")
                    continue
                for i, s in enumerate(series):
                    many = len(series) > 1
                    images.append(ImageEntry(str(p), f"{p.stem} · {s}" if many else p.stem, cond.name, str(p.parent),
                                             i, f"{p}::{i}" if many else str(p)))
    return images, remarks
