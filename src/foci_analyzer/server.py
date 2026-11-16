"""The window: a small web server on this computer only, and the pages in `static/`.

Requests under /api/ need the token that is written into the page when it is served, and a Host header of
this computer, so that other web pages cannot use the server.
"""

from __future__ import annotations

import json
import mimetypes
import os
import secrets
import subprocess
import sys
import threading
import urllib.parse
from collections import OrderedDict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from pathlib import Path

import numpy as np

from . import __version__, batch, config, experiment, models, pipeline, stats
from .imagestore import ImageStore
from .stages import StageCache

_PICK = """
import sys, tkinter as tk
from tkinter import filedialog
root = tk.Tk(); root.withdraw(); root.attributes("-topmost", True)
if sys.argv[1] == "files":
    out = filedialog.askopenfilenames(title="Choose images", filetypes=[("Images", "*.tif *.tiff *.czi *.lif"), ("All files", "*.*")])
    sys.stdout.write("\\n".join(out))
else:
    sys.stdout.write(filedialog.askdirectory(title="Choose a folder") or "")
"""


def home() -> Path:
    """Where the window keeps its last settings and experiment."""
    return Path(os.environ.get("FOCI_ANALYZER_HOME") or Path.home() / ".foci_analyzer")


class ApiError(Exception):
    def __init__(self, message: str, status: int = 400):
        super().__init__(message)
        self.status = status


def _wrap16(labels: np.ndarray) -> np.ndarray:
    """Labels as 16-bit for the viewer (it only needs neighbours to differ)."""
    if labels.max(initial=0) < 65536:
        return labels.astype("<u2")
    return np.where(labels > 0, (labels - 1) % 65535 + 1, 0).astype("<u2")


def _columns(df) -> dict:
    """A table as {column: values} for JSON (NaN becomes null)."""
    d = json.loads(df.to_json(orient="split", index=False))
    return {c: [row[j] for row in d["data"]] for j, c in enumerate(d["columns"])}


class App:
    """What the window works on: the settings, the experiment, the opened images and the results."""

    def __init__(self):
        self.token = secrets.token_urlsafe(16)
        self.lock = threading.Lock()
        self.compute_lock = threading.Lock()  # one analysis at a time
        self.store = ImageStore()
        self.cache = StageCache(keep=40)
        self.settings = config.Settings()
        self.conditions: list[experiment.Condition] = []
        self.images: list[experiment.ImageEntry] = []
        self.remarks: list[str] = []
        self.output = ""
        self.results = batch.Results()
        self.progress = batch.Progress()
        self._planes: OrderedDict[tuple, np.ndarray] = OrderedDict()
        self._load()

    # ---- state on disk
    def _load(self) -> None:
        try:
            self.settings = config.load(home() / "last_settings.toml")
        except (OSError, ValueError, TypeError):
            pass  # no or unusable file: start from the defaults
        try:
            saved = json.loads((home() / "last_experiment.json").read_text(encoding="utf-8"))
            if isinstance(saved, list):
                saved = {"conditions": saved}
            self.set_conditions(saved.get("conditions", []), saved.get("output", ""), save=False)
            self.results.load(self.results_folder())
        except (OSError, ValueError, ApiError):
            pass

    def set_settings(self, data: dict) -> dict:
        try:
            settings = config.from_dict(data)
        except (ValueError, TypeError, AttributeError) as e:
            raise ApiError(str(e)) from e
        with self.lock:
            self.settings = settings
        config.save(settings, home() / "last_settings.toml")
        return config.to_dict(settings)

    def set_conditions(self, items: list[dict], output: str = "", save: bool = True) -> dict:
        try:
            conditions = experiment.conditions_from_dicts(items)
        except (ValueError, TypeError, AttributeError) as e:
            raise ApiError(str(e)) from e
        images, remarks = experiment.scan(conditions)
        with self.lock:
            self.conditions, self.images, self.remarks = conditions, images, remarks
            self.output = str(output or "").strip().strip('"')
        if save:
            home().mkdir(parents=True, exist_ok=True)
            (home() / "last_experiment.json").write_text(
                json.dumps({"conditions": self._conditions(), "output": self.output}), encoding="utf-8")
        return self.experiment_state()

    def results_folder(self) -> str:
        """The chosen results folder, or `Foci_Analyzer_results` next to the first image."""
        if self.output:
            return self.output
        return str(Path(self.images[0].folder) / "Foci_Analyzer_results") if self.images else ""

    def _conditions(self) -> list[dict]:
        return [vars(c) for c in self.conditions]

    def experiment_state(self) -> dict:
        return {"conditions": self._conditions(), "images": [i.to_dict() for i in self.images],
                "remarks": self.remarks, "output": self.output, "results_folder": self.results_folder()}

    def state(self) -> dict:
        return {"version": __version__, "schema": config.schema(), "settings": config.to_dict(self.settings),
                "models": models.available(), "batch": self.progress.to_dict(), **self.experiment_state()}

    # ---- images
    def entry(self, image_id: str) -> experiment.ImageEntry:
        for i in self.images:  # only images of the experiment can be read
            if i.id == image_id:
                return i
        raise ApiError("this image is not in the experiment", 404)

    def image(self, image_id: str):
        e = self.entry(image_id)
        try:
            return self.store.get(e.id, e.path, e.series)
        except (OSError, ValueError) as err:
            raise ApiError(f"cannot read {Path(e.path).name}: {err}") from err

    def image_info(self, image_id: str) -> dict:
        im = self.image(image_id)
        return {**im.info.to_dict(), "ranges": [im.display_range(c) for c in range(im.info.channels)],
                "scale": im.scale}

    def _plane_args(self, q: dict, info):
        try:
            t, level = int(q.get("t", 0)), int(q.get("level", 0))
            z = None if q.get("z", "mip") == "mip" else int(q["z"])
        except (KeyError, ValueError) as e:
            raise ApiError("t, z and level must be numbers (z may be 'mip')") from e
        if not (0 <= t < info.frames and 0 <= level <= 8 and (z is None or 0 <= z < info.slices)):
            raise ApiError("outside the image")
        return t, z, level

    @staticmethod
    def _crop(plane: np.ndarray, q: dict) -> tuple[bytes, int, int]:
        h, w = plane.shape
        try:
            x0, y0 = int(q.get("x", 0)), int(q.get("y", 0))
            x1, y1 = x0 + int(q.get("w", w)), y0 + int(q.get("h", h))
        except ValueError as e:
            raise ApiError("x, y, w and h must be numbers") from e
        x0, y0, x1, y1 = max(0, x0), max(0, y0), min(w, x1), min(h, y1)
        if x1 <= x0 or y1 <= y0:
            raise ApiError("outside the image")
        return np.ascontiguousarray(plane[y0:y1, x0:x1]).astype("<u2", copy=False).tobytes(), x1 - x0, y1 - y0

    def region(self, q: dict) -> tuple[bytes, int, int]:
        im = self.image(q.get("path", ""))
        t, z, level = self._plane_args(q, im.info)
        try:
            c = int(q["c"])
        except (KeyError, ValueError) as e:
            raise ApiError("c must be a number") from e
        if not 0 <= c < im.info.channels:
            raise ApiError("outside the image")
        return self._crop(im.plane(t, c, z, level), q)

    # ---- analysis of the image in the viewer
    def frame(self, image_id: str, t: int, only_candidates: bool = False) -> pipeline.FrameResult:
        im = self.image(image_id)
        if not 0 <= t < im.info.frames:
            raise ApiError("outside the image")
        with self.lock:
            settings = self.settings
        try:
            with self.compute_lock:
                return pipeline.analyze_frame(self.cache, f"{image_id}|{t}", im.data[t], im.info, settings,
                                              only_candidates)
        except ValueError as e:
            raise ApiError(str(e)) from e

    def preview(self, body: dict) -> dict:
        if "settings" in body:
            self.set_settings(body["settings"])
        r = self.frame(body.get("path", ""), int(body.get("t", 0)))
        return {"key": r.key, "cells": int(r.labels.max()), "channels": [c.summary() for c in r.channels],
                "remarks": r.remarks, "pixel_size": r.cal.pixel_size, "binning": r.binning, "border": r.border,
                "table": _columns(r.cells)}

    def candidates(self, q: dict) -> bytes:
        """Per candidate: x, y, z (image pixels), height, cell, and 1 when it is a focus (float32)."""
        r = self.frame(q.get("path", ""), int(q.get("t", 0)))
        ch = next((c for c in r.channels if str(c.channel) == q.get("ch")), None)
        if ch is None:
            raise ApiError("no such foci channel")
        out = np.zeros((len(ch.height), 6), "<f4")
        off = r.border + (r.binning - 1) / 2
        out[:, 0], out[:, 1] = ch.pos[:, 2] * r.binning + off, ch.pos[:, 1] * r.binning + off
        out[:, 2], out[:, 3], out[:, 4] = ch.pos[:, 0], ch.height, ch.cell
        out[ch.kept, 5] = 1
        if ch.labels is not None and ch.labels.shape[0] == 1:
            out[:, 2] = -1  # found in a projection: shown in every slice
        return out.tobytes()

    def overlay(self, q: dict) -> tuple[bytes, int, int]:
        """A region of the label image of the cells or of the foci of one channel, in image pixels."""
        im = self.image(q.get("path", ""))
        t, z, level = self._plane_args(q, im.info)
        r = self.frame(q.get("path", ""), t)
        kind = q.get("kind", "cells")
        key = (r.key, kind, q.get("ch"), z)
        with self.lock:
            plane = self._planes.get(key)
        if plane is None:
            if kind == "cells":
                lab = r.labels
            else:
                ch = next((c for c in r.channels if str(c.channel) == q.get("ch")), None)
                if ch is None:
                    raise ApiError("no such foci channel")
                lab = ch.labels
            small = lab[0] if lab.shape[0] == 1 else (lab.max(0) if z is None else lab[z])
            if r.binning > 1:
                small = np.repeat(np.repeat(small, r.binning, 0), r.binning, 1)
            plane = np.zeros((im.info.height, im.info.width), "<u2")
            plane[r.border:r.border + small.shape[0], r.border:r.border + small.shape[1]] = _wrap16(small)
            with self.lock:
                self._planes[key] = plane
                while len(self._planes) > 12:
                    self._planes.popitem(last=False)
        f = 1 << level
        return self._crop(plane[::f, ::f], q)

    # ---- batch and results
    def start_batch(self, body: dict) -> dict:
        if "settings" in body:
            self.set_settings(body["settings"])
        if self.progress.running:
            raise ApiError("an analysis is already running")
        if not self.images:
            raise ApiError("there are no images; add them on the Experiment page")
        wanted = set(body.get("images") or [])
        images = [i for i in self.images if not wanted or i.id in wanted]
        self.progress.running = True
        threading.Thread(target=batch.run, daemon=True, args=(
            images, self.settings, self.results_folder(), self.results, self.progress, self.compute_lock)).start()
        return self.progress.to_dict()

    def results_json(self) -> bytes:
        r = self.results
        meta = dict(r.meta)
        meta["current"] = bool(r.meta) and r.meta.get("settings") == config.to_dict(self.settings)
        meta["conditions"] = [c.name for c in self.conditions]
        return (f'{{"cells":{r.cells.to_json(orient="split", index=False) if len(r.cells) else "null"},'
                f'"foci":{r.foci.to_json(orient="split", index=False) if len(r.foci) else "null"},'
                f'"meta":{json.dumps(meta)}}}').encode()

    def compare(self, body: dict) -> dict:
        df = self.results.foci if body.get("table") == "foci" else self.results.cells
        if not len(df):
            raise ApiError("there are no results yet")
        if body.get("table") == "foci" and "folder" not in df:
            raise ApiError("the foci table has no folder column")
        keep = body.get("conditions")
        if keep:
            df = df[df["condition"].isin(keep)]
        if body.get("table") == "foci" and body.get("channel"):
            df = df[df["channel"] == body["channel"]]
        try:
            return stats.compare(df, str(body.get("measure", "")), str(body.get("unit", "cell")))
        except ValueError as e:
            raise ApiError(str(e)) from e

    def open_results_folder(self) -> dict:
        folder = self.results.meta.get("folder") or self.results_folder()
        if not folder or not Path(folder).is_dir():
            raise ApiError("there is no results folder yet")
        os.startfile(folder)
        return {"folder": folder}


def pick(kind: str) -> list[str]:
    """Folders or files chosen in a normal dialog (a browser cannot give the path of a folder)."""
    try:
        out = subprocess.run([sys.executable, "-c", _PICK, kind], capture_output=True, text=True, timeout=600,
                             check=False)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise ApiError(f"the dialog could not be opened: {e}", 500) from e
    if out.returncode != 0:
        raise ApiError("the dialog could not be opened; type or paste the folder instead", 500)
    return [line for line in out.stdout.splitlines() if line.strip()]


def make_handler(app: App):
    static = resources.files("foci_analyzer") / "static"

    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, *args):  # quiet
            pass

        # ---- answers
        def _send(self, status: int, body: bytes, ctype: str, headers: dict | None = None):
            self.send_response(status)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            for k, v in (headers or {}).items():
                self.send_header(k, str(v))
            self.end_headers()
            self.wfile.write(body)

        def _json(self, data, status: int = 200):
            self._send(status, json.dumps(data).encode(), "application/json")

        def _allowed(self) -> bool:
            host = (self.headers.get("Host") or "").split(":")[0]
            return host in ("127.0.0.1", "localhost") and secrets.compare_digest(
                self.headers.get("X-Token") or "", app.token)

        def _api(self, call):
            if not self._allowed():
                return self._json({"error": "not allowed"}, 403)
            try:
                call()
            except ApiError as e:
                self._json({"error": str(e)}, e.status)
            except Exception as e:  # noqa: BLE001  a bug: report it in the window, not as a broken connection
                self._json({"error": f"{type(e).__name__}: {e}"}, 500)

        def _binary(self, result):
            if isinstance(result, tuple):
                body, w, h = result
                return self._send(200, body, "application/octet-stream", {"X-Width": w, "X-Height": h})
            self._send(200, result, "application/octet-stream")

        # ---- requests
        def do_GET(self):
            url = urllib.parse.urlsplit(self.path)
            q = dict(urllib.parse.parse_qsl(url.query))
            if url.path == "/":
                page = (static / "app.html").read_text(encoding="utf-8").replace("{{TOKEN}}", app.token)
                return self._send(200, page.encode(), "text/html; charset=utf-8")
            if url.path == "/static/plotly.min.js":
                import plotly
                body = (Path(plotly.__file__).parent / "package_data" / "plotly.min.js").read_bytes()
                return self._send(200, body, "text/javascript", {"Cache-Control": "max-age=86400"})
            if url.path.startswith("/static/"):
                name = url.path[len("/static/"):]
                if "/" in name or "\\" in name or name.startswith(".") or not (static / name).is_file():
                    return self._send(404, b"not found", "text/plain")
                ctype = mimetypes.guess_type(name)[0] or "application/octet-stream"
                if name.endswith(".js"):
                    ctype = "text/javascript"
                return self._send(200, (static / name).read_bytes(), ctype)
            gets = {
                "/api/state": lambda: self._json(app.state()),
                "/api/image/info": lambda: self._json(app.image_info(q.get("path", ""))),
                "/api/image/region": lambda: self._binary(app.region(q)),
                "/api/overlay": lambda: self._binary(app.overlay(q)),
                "/api/candidates": lambda: self._binary(app.candidates(q)),
                "/api/batch": lambda: self._json(app.progress.to_dict()),
                "/api/results": lambda: self._send(200, app.results_json(), "application/json"),
            }
            if url.path in gets:
                return self._api(gets[url.path])
            self._send(404, b"not found", "text/plain")

        def do_POST(self):
            path = urllib.parse.urlsplit(self.path).path
            try:
                body = json.loads(self.rfile.read(int(self.headers.get("Content-Length") or 0)) or b"{}")
            except ValueError:
                return self._json({"error": "the request is not JSON"}, 400)

            def stop():
                app.progress.stop = True
                return app.progress.to_dict()

            posts = {
                "/api/settings": lambda: app.set_settings(body.get("settings", {})),
                "/api/experiment": lambda: app.set_conditions(body.get("conditions", []), body.get("output", "")),
                "/api/pick": lambda: {"paths": pick("files" if body.get("kind") == "files" else "folder")},
                "/api/preview": lambda: app.preview(body),
                "/api/batch/start": lambda: app.start_batch(body),
                "/api/batch/stop": stop,
                "/api/stats": lambda: app.compare(body),
                "/api/results/open": app.open_results_folder,
            }
            if path in posts:
                return self._api(lambda: self._json(posts[path]()))
            self._send(404, b"not found", "text/plain")

    return Handler


def make_server(port: int = 8770, app: App | None = None) -> tuple[ThreadingHTTPServer, App]:
    app = app or App()
    server = ThreadingHTTPServer(("127.0.0.1", port), make_handler(app))
    server.daemon_threads = True
    return server, app
