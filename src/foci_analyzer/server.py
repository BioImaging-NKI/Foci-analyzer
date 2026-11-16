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
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from pathlib import Path

from . import __version__, config, experiment
from .imagestore import ImageStore

_PICK_FOLDER = """
import sys, tkinter as tk
from tkinter import filedialog
root = tk.Tk(); root.withdraw(); root.attributes("-topmost", True)
sys.stdout.write(filedialog.askdirectory(title="Choose a folder with images") or "")
"""


def home() -> Path:
    """Where the window keeps its last settings and experiment."""
    return Path(os.environ.get("FOCI_ANALYZER_HOME") or Path.home() / ".foci_analyzer")


class ApiError(Exception):
    def __init__(self, message: str, status: int = 400):
        super().__init__(message)
        self.status = status


class App:
    """What the window works on: the settings, the experiment and the opened images."""

    def __init__(self):
        self.token = secrets.token_urlsafe(16)
        self.lock = threading.Lock()
        self.store = ImageStore()
        self.settings = config.Settings()
        self.conditions: list[experiment.Condition] = []
        self.images: list[experiment.ImageEntry] = []
        self.remarks: list[str] = []
        self._load()

    # ---- state on disk
    def _load(self) -> None:
        try:
            self.settings = config.load(home() / "last_settings.toml")
        except (OSError, ValueError, TypeError):
            pass  # no or unusable file: start from the defaults
        try:
            items = json.loads((home() / "last_experiment.json").read_text(encoding="utf-8"))
            self.set_conditions(items, save=False)
        except (OSError, ValueError):
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

    def set_conditions(self, items: list[dict], save: bool = True) -> dict:
        try:
            conditions = experiment.conditions_from_dicts(items)
        except (ValueError, TypeError, AttributeError) as e:
            raise ApiError(str(e)) from e
        images, remarks = experiment.scan(conditions)
        with self.lock:
            self.conditions, self.images, self.remarks = conditions, images, remarks
        if save:
            home().mkdir(parents=True, exist_ok=True)
            (home() / "last_experiment.json").write_text(json.dumps(self._conditions()), encoding="utf-8")
        return self.experiment_state()

    def _conditions(self) -> list[dict]:
        return [vars(c) for c in self.conditions]

    def experiment_state(self) -> dict:
        return {"conditions": self._conditions(), "images": [i.to_dict() for i in self.images],
                "remarks": self.remarks}

    def state(self) -> dict:
        return {"version": __version__, "schema": config.schema(), "settings": config.to_dict(self.settings),
                **self.experiment_state()}

    # ---- images
    def image(self, path: str):
        if path not in {i.path for i in self.images}:  # only images of the experiment can be read
            raise ApiError("this image is not in the experiment", 404)
        try:
            return self.store.get(path)
        except (OSError, ValueError) as e:
            raise ApiError(f"cannot read {Path(path).name}: {e}") from e

    def image_info(self, path: str) -> dict:
        im = self.image(path)
        return {**im.info.to_dict(), "ranges": [im.display_range(c) for c in range(im.info.channels)]}

    def region(self, q: dict) -> tuple[bytes, int, int]:
        im = self.image(q.get("path", ""))
        try:
            t, c, level = int(q.get("t", 0)), int(q["c"]), int(q.get("level", 0))
            z = None if q.get("z", "mip") == "mip" else int(q["z"])
        except (KeyError, ValueError) as e:
            raise ApiError("region: t, c, z and level must be numbers (z may be 'mip')") from e
        info = im.info
        if not (0 <= t < info.frames and 0 <= c < info.channels and 0 <= level <= 8
                and (z is None or 0 <= z < info.slices)):
            raise ApiError("region: outside the image")
        plane = im.plane(t, c, z, level)
        h, w = plane.shape
        try:
            x0, y0 = int(q.get("x", 0)), int(q.get("y", 0))
            x1, y1 = x0 + int(q.get("w", w)), y0 + int(q.get("h", h))
        except ValueError as e:
            raise ApiError("region: x, y, w and h must be numbers") from e
        x0, y0, x1, y1 = max(0, x0), max(0, y0), min(w, x1), min(h, y1)
        if x1 <= x0 or y1 <= y0:
            raise ApiError("region: outside the image")
        part = plane[y0:y1, x0:x1]
        return part.astype("<u2", copy=False).tobytes(), x1 - x0, y1 - y0


def pick_folder() -> str:
    """A folder chosen in a normal folder dialog (a browser cannot give the path of a folder)."""
    try:
        out = subprocess.run([sys.executable, "-c", _PICK_FOLDER], capture_output=True, text=True, timeout=600, check=False)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise ApiError(f"the folder dialog could not be opened: {e}", 500) from e
    if out.returncode != 0:
        raise ApiError("the folder dialog could not be opened; type or paste the folder instead", 500)
    return out.stdout.strip()


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

        # ---- requests
        def do_GET(self):
            url = urllib.parse.urlsplit(self.path)
            q = dict(urllib.parse.parse_qsl(url.query))
            if url.path == "/":
                page = (static / "app.html").read_text(encoding="utf-8").replace("{{TOKEN}}", app.token)
                return self._send(200, page.encode(), "text/html; charset=utf-8")
            if url.path.startswith("/static/"):
                name = url.path[len("/static/"):]
                if "/" in name or "\\" in name or name.startswith(".") or not (static / name).is_file():
                    return self._send(404, b"not found", "text/plain")
                ctype = mimetypes.guess_type(name)[0] or "application/octet-stream"
                if name.endswith(".js"):
                    ctype = "text/javascript"
                return self._send(200, (static / name).read_bytes(), ctype)
            if url.path == "/api/state":
                return self._api(lambda: self._json(app.state()))
            if url.path == "/api/image/info":
                return self._api(lambda: self._json(app.image_info(q.get("path", ""))))
            if url.path == "/api/image/region":
                def region():
                    body, w, h = app.region(q)
                    self._send(200, body, "application/octet-stream", {"X-Width": w, "X-Height": h})
                return self._api(region)
            self._send(404, b"not found", "text/plain")

        def do_POST(self):
            path = urllib.parse.urlsplit(self.path).path
            try:
                body = json.loads(self.rfile.read(int(self.headers.get("Content-Length") or 0)) or b"{}")
            except ValueError:
                return self._json({"error": "the request is not JSON"}, 400)
            if path == "/api/settings":
                return self._api(lambda: self._json(app.set_settings(body.get("settings", {}))))
            if path == "/api/experiment":
                return self._api(lambda: self._json(app.set_conditions(body.get("conditions", []))))
            if path == "/api/pick-folder":
                return self._api(lambda: self._json({"folder": pick_folder()}))
            self._send(404, b"not found", "text/plain")

    return Handler


def make_server(port: int = 8770, app: App | None = None) -> tuple[ThreadingHTTPServer, App]:
    app = app or App()
    server = ThreadingHTTPServer(("127.0.0.1", port), make_handler(app))
    server.daemon_threads = True
    return server, app
