"""Command line: starts the window in the browser.

The port is taken and the browser opened before the program's packages are loaded (numpy, scipy, ... take
seconds on a slow computer): until then the browser shows a starting page with the logo, which reloads
itself into the window as soon as the program answers.
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
import threading
import time
import traceback
import urllib.request
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from pathlib import Path

from . import __version__

STARTING = """<!doctype html><html><head><meta charset="utf-8"><title>Foci Analyzer (starting)</title>
<link rel="icon" href="/static/favicon.png"><style>
body{margin:0;height:100vh;display:flex;align-items:center;justify-content:center;font:16px system-ui,sans-serif;
background:#f4f5f7;color:#333}@media(prefers-color-scheme:dark){body{background:#1d1f23;color:#ddd}}
div.box{text-align:center}.logo{position:relative;width:168px;height:168px;margin:0 auto}
.logo img{position:absolute;left:20px;top:20px;width:128px;height:128px;animation:p 1.6s ease-in-out infinite}
.logo svg{position:absolute;left:0;top:0;transform:rotate(-90deg)}
.logo circle{fill:none;stroke-width:5}.logo .track{stroke:rgba(127,127,127,.18)}
.logo .ring{stroke:#1f7a4d;stroke-linecap:round;transition:stroke-dashoffset .35s linear}
@media(prefers-color-scheme:dark){.logo .ring{stroke:#4fb783}}
@keyframes p{50%{opacity:.45}}p{margin:14px 0 0}small{opacity:.7}</style></head>
<body data-starting><div class="box"><div class="logo"><svg width="168" height="168"><circle class="track" cx="84" cy="84" r="80"/>
<circle class="ring" id="ring" cx="84" cy="84" r="80" stroke-dasharray="502.7" stroke-dashoffset="502.7"/></svg>
<img src="/static/logo.png" alt=""></div><p>Starting Foci Analyzer&hellip;</p>
<small id="step">The window opens here in a moment.</small></div><script>
// reloads once the answer is no longer this page (the window, or the error); the page must not contain the
// text it looks for, else it reloads itself all the time
(function poll(){fetch("/",{cache:"no-store"}).then(r=>r.text()).then(t=>{
if(t.indexOf("data-starting")<0)location.reload();else setTimeout(poll,400);})
.catch(()=>setTimeout(poll,1000));})();
// the ring around the logo: how far the start is (the steps, timed as at the last start)
(function step(){fetch("/starting.json",{cache:"no-store"}).then(r=>r.json()).then(s=>{
document.getElementById("ring").style.strokeDashoffset=String(502.7*(1-s.fraction));
if(s.text)document.getElementById("step").textContent=s.text+"…";setTimeout(step,150);})
.catch(()=>setTimeout(step,500));})();</script></body></html>"""


class Startup:
    """The steps of a start, timed: the starting page shows how far it is, by the times of the last start."""

    def __init__(self, folder: Path):
        self.file = folder / "startup.json"
        try:
            self.last = [float(x) for x in json.loads(self.file.read_text(encoding="utf-8"))["seconds"]]
        except (OSError, ValueError, KeyError, TypeError):
            self.last = []
        self.seconds: list[float] = []
        self.text, self.t = "", time.perf_counter()

    def next(self, text: str) -> None:
        now = time.perf_counter()
        if self.text:
            self.seconds.append(now - self.t)
        self.text, self.t = text, now

    def fraction(self) -> float:
        n = len(STEPS) + 2  # the packages, then the program itself, then the settings
        plan = self.last if len(self.last) == n else [1.0] * n  # the first start: every step alike
        i = len(self.seconds)
        done = sum(plan[:i]) + (min(time.perf_counter() - self.t, 0.95 * plan[i]) if i < n else 0.0)
        return max(0.02, min(0.98, done / sum(plan)))

    def finish(self) -> None:
        self.next("")
        try:
            self.file.parent.mkdir(parents=True, exist_ok=True)
            self.file.write_text(json.dumps({"seconds": [round(x, 3) for x in self.seconds]}), encoding="utf-8")
        except OSError:
            pass


# the packages, loaded one by one so that the starting page can say what it is doing (the slow part of a start)
STEPS = [("numpy", "Loading the numerical libraries"), ("scipy.ndimage", "Loading the image filters"),
         ("skimage.segmentation", "Loading the image analysis"), ("pandas", "Loading the tables"),
         ("scipy.stats", "Loading the statistics"), ("tifffile", "Loading the image readers")]


class LocalServer(ThreadingHTTPServer):
    # Windows lets a second program take a port that is in use when SO_REUSEADDR is set: two Foci Analyzers
    # (say an installed one and a development one) then both answer on 8770, at random.
    allow_reuse_address = sys.platform != "win32"
    daemon_threads = True


class StartingHandler(BaseHTTPRequestHandler):
    """Answers while the program loads: the starting page, and the logo it shows."""

    error = ""  # set when the program could not start: shown instead of the starting page
    startup: Startup | None = None  # how far the start is (for the ring on the starting page)
    seen = False  # a browser asked for something (after an update: the tab that started it is still open)

    def log_message(self, *args):
        pass

    def do_GET(self):
        StartingHandler.seen = True
        name = self.path.split("?")[0].rsplit("/", 1)[-1]
        if self.path.split("?")[0] == "/starting.json":
            st = StartingHandler.startup
            body = json.dumps({"fraction": round(st.fraction(), 3) if st else 0.02,
                               "text": st.text if st else ""}).encode()
            kind = "application/json"
        elif self.path.startswith("/static/") and name in ("logo.png", "favicon.png"):
            body, kind = (resources.files("foci_analyzer") / "static" / name).read_bytes(), "image/png"
        elif self.path.startswith("/api/"):
            self.send_response(503)
            self.end_headers()
            return
        elif self.error:
            text = self.error.replace("&", "&amp;").replace("<", "&lt;")
            body = (f"<!doctype html><title>Foci Analyzer (error)</title><h3>Foci Analyzer could not start</h3>"
                    f"<pre>{text}</pre>").encode()
            kind = "text/html; charset=utf-8"
        else:
            body, kind = STARTING.encode(), "text/html; charset=utf-8"
        self.send_response(200)
        self.send_header("Content-Type", kind)
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def _already_running(url: str) -> bool:
    try:
        with urllib.request.urlopen(url, timeout=2) as r:
            return b"<title>Foci Analyzer" in r.read(600)
    except OSError:
        return False


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="foci-analyzer", description="Foci Analyzer: opens the window in the browser.")
    p.add_argument("--port", type=int, default=8770)
    p.add_argument("--no-browser", action="store_true", help="start the server without opening a browser tab")
    p.add_argument("--home", help="folder for the last settings and experiment (default: ~/.foci_analyzer); "
                   "a separate one gives a test window that leaves your own state alone")
    p.add_argument("--after-update", action="store_true",
                   help="started by an update: open the browser only if no tab is waiting for the window")
    p.add_argument("--version", action="version", version=__version__)
    args = p.parse_args(argv)
    if args.home:
        os.environ["FOCI_ANALYZER_HOME"] = args.home

    url = f"http://127.0.0.1:{args.port}/"
    try:
        httpd = LocalServer(("127.0.0.1", args.port), StartingHandler)
    except OSError as e:
        if _already_running(url):  # started twice: show the window that is there
            if not args.no_browser:
                webbrowser.open(url)
            return 0
        print(f"Port {args.port} cannot be used ({e}). Try another one with --port.")
        return 1
    serving = threading.Thread(target=httpd.serve_forever, daemon=True)
    serving.start()
    if not args.no_browser and not args.after_update:
        webbrowser.open(url)

    startup = StartingHandler.startup = Startup(Path(os.environ.get("FOCI_ANALYZER_HOME") or Path.home() / ".foci_analyzer"))
    try:
        for module, text in STEPS:
            startup.next(text)
            importlib.import_module(module)
        startup.next("Loading the analysis and the window")
        from . import models
        from .server import App, make_handler

        startup.next("Reading your settings and experiment")
        app = App()
        app.warm = True
        app.port = httpd.server_address[1]
        app.stop_server = httpd.shutdown
        app.warm_up()
        app.check_updates_soon()
    except Exception:  # noqa: BLE001 - shown in the browser: there may be no console
        StartingHandler.error = traceback.format_exc()
        print(StartingHandler.error)
        serving.join(60)  # long enough for the browser to show it
        httpd.shutdown()
        return 1
    startup.finish()
    httpd.RequestHandlerClass = make_handler(app)  # read per request: from now on the window answers
    print(f"Foci Analyzer {__version__} at {url}  (Ctrl+C to stop)")
    if args.after_update and not args.no_browser:
        time.sleep(3)  # the tab that started the update asks every second
        if not StartingHandler.seen:
            webbrowser.open(url)
    try:
        while serving.is_alive():  # /api/quit stops serve_forever
            serving.join(0.5)
    except KeyboardInterrupt:
        httpd.shutdown()
    finally:
        httpd.server_close()
        models.stop_all()
    return 0
