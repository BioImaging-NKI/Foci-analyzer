"""Analyze all images of the experiment and keep the tables (plan.md, section 7).

The results folder holds `cells.csv` (one row per cell and frame), `foci.csv` (one row per focus),
`settings.toml` and `run.json`. Files are written under a temporary name and then renamed, so an interrupted
run never leaves a half-written table.
"""

from __future__ import annotations

import json
import os
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import __version__, config, pipeline, tracking
from .experiment import ImageEntry
from .io import read_image
from .stages import StageCache


@dataclass
class Progress:
    running: bool = False
    done: int = 0
    total: int = 0
    current: str = ""
    errors: list[str] = field(default_factory=list)
    stop: bool = False
    finished: float = 0.0

    def to_dict(self) -> dict:
        return {"running": self.running, "done": self.done, "total": self.total, "current": self.current,
                "errors": self.errors, "finished": self.finished}


def _replace(path: Path, write) -> None:
    tmp = path.with_name(path.name + ".tmp")
    write(tmp)
    os.replace(tmp, path)


def analyze_image(entry: ImageEntry, s: config.Settings, compute_lock: threading.Lock | None = None):
    """The tables of one image (all frames), with the columns that say where a row comes from."""
    data, info = read_image(entry.path, entry.series)
    cache = StageCache(keep=24)
    cells, foci, labels, frames = [], [], [], []
    ch = s.image.cell_channel or s.image.nuclei_channel
    remark = ""
    for t in range(info.frames):
        if compute_lock:
            with compute_lock:
                r = pipeline.analyze_frame(cache, f"{entry.path}|{t}", data[t], info, s)
        else:
            r = pipeline.analyze_frame(cache, f"{entry.path}|{t}", data[t], info, s)
        for df, out in ((r.cells, cells), (r.foci, foci)):
            df = df.copy()
            df.insert(0, "frame", t + 1)
            out.append(df)
        labels.append(r.labels)
        if ch and s.timelapse.tracking == "trackastra":
            frames.append(data[t, ch - 1].max(0))
    cells, foci = pd.concat(cells, ignore_index=True), pd.concat(foci, ignore_index=True)
    if info.frames > 1 and s.timelapse.tracking != "off":
        remark = tracking.add_tracks(cells, foci, labels, s.timelapse, frames or None)
    for df in (cells, foci):
        for i, (col, value) in enumerate((("condition", entry.condition), ("image", entry.name),
                                          ("folder", entry.folder))):
            df.insert(i, col, value)
    cells["path"] = entry.id
    meta = {"unit": info.unit, "pixel_size": r.cal.pixel_size, "frames": info.frames, "slices": info.slices,
            "frame_interval": info.frame_interval, "remark": remark}
    return cells, foci, meta


class Results:
    """The tables of the last run, in memory and in the results folder."""

    def __init__(self):
        self.cells, self.foci, self.meta = pd.DataFrame(), pd.DataFrame(), {}

    def load(self, folder: str | Path) -> bool:
        folder = Path(folder)
        try:
            self.meta = json.loads((folder / "run.json").read_text(encoding="utf-8"))
            self.cells = pd.read_csv(folder / "cells.csv", keep_default_na=True)
            self.cells["check"] = self.cells["check"].fillna("") if "check" in self.cells else ""
            try:
                self.foci = pd.read_csv(folder / "foci.csv")
            except pd.errors.EmptyDataError:
                self.foci = pd.DataFrame()
        except (OSError, ValueError, KeyError):
            return False
        return True

    def save(self, folder: str | Path, s: config.Settings) -> None:
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        _replace(folder / "cells.csv", lambda p: self.cells.to_csv(p, index=False))
        _replace(folder / "foci.csv", lambda p: self.foci.to_csv(p, index=False))
        _replace(folder / "settings.toml", lambda p: config.save(s, p))
        _replace(folder / "run.json", lambda p: p.write_text(json.dumps(self.meta, indent=1), encoding="utf-8"))


def run(images: list[ImageEntry], s: config.Settings, folder: str, results: Results, progress: Progress,
        compute_lock: threading.Lock) -> None:
    """Analyze the images one after another (meant for a thread); an image that fails is reported and skipped."""
    progress.running, progress.done, progress.total, progress.errors, progress.stop = True, 0, len(images), [], False
    cells, foci, per_image = [], [], {}
    try:
        # the networks share one GPU, so images are only analyzed side by side without them
        networks = s.segmentation.method in ("cellpose", "stardist", "instanseg") or any(
            f.method == "spotiflow" for f in s.foci)
        workers = 1 if networks else min(len(images), max(1, (os.cpu_count() or 2) // 2), 6)
        done: dict[int, tuple] = {}

        def finished(i, entry, get):
            try:
                done[i] = get()
            except Exception as e:  # noqa: BLE001  one unreadable image must not stop the run
                progress.errors.append(f"{entry.name}: {e}")
            progress.done += 1

        if workers <= 1:
            for i, entry in enumerate(images):
                if progress.stop:
                    break
                progress.current = entry.name
                finished(i, entry, lambda entry=entry: analyze_image(entry, s, compute_lock))
        else:
            progress.current = f"{workers} images at a time"
            with ProcessPoolExecutor(workers) as pool:
                futures = {pool.submit(analyze_image, entry, s): (i, entry) for i, entry in enumerate(images)}
                for fut in as_completed(futures):
                    finished(*futures[fut], fut.result)
                    if progress.stop:
                        pool.shutdown(wait=False, cancel_futures=True)
                        break
        if progress.stop:
            progress.errors.append("stopped before all images were analyzed")
        for i in sorted(done):
            cells.append(done[i][0])
            foci.append(done[i][1])
            per_image[images[i].id] = done[i][2]
        if cells:
            results.cells = pd.concat(cells, ignore_index=True)
            results.foci = pd.concat(foci, ignore_index=True)
            units = {m["unit"] for m in per_image.values()}
            results.meta = {"version": __version__, "when": time.strftime("%Y-%m-%d %H:%M"), "folder": str(folder),
                            "settings": config.to_dict(s), "unit": units.pop() if len(units) == 1 else "mixed",
                            "images": per_image, "errors": progress.errors}
            try:
                results.save(folder, s)
            except OSError as e:
                progress.errors.append(f"the results could not be saved in {folder}: {e}")
    finally:
        progress.running, progress.current, progress.finished = False, "", time.time()
