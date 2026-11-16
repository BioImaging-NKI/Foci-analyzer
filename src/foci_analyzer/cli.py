"""Command line: starts the window in the browser."""

from __future__ import annotations

import argparse
import os
import urllib.request
import webbrowser

from . import __version__, models
from .server import make_server


def _already_running(url: str) -> bool:
    try:
        with urllib.request.urlopen(url, timeout=2) as r:
            return b"<title>Foci Analyzer</title>" in r.read(600)
    except OSError:
        return False


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="foci-analyzer", description="Foci Analyzer: opens the window in the browser.")
    p.add_argument("--port", type=int, default=8770)
    p.add_argument("--no-browser", action="store_true", help="start the server without opening a browser tab")
    p.add_argument("--home", help="folder for the last settings and experiment (default: ~/.foci_analyzer); "
                   "a separate one gives a test window that leaves your own state alone")
    p.add_argument("--version", action="version", version=__version__)
    args = p.parse_args(argv)
    if args.home:
        os.environ["FOCI_ANALYZER_HOME"] = args.home

    url = f"http://127.0.0.1:{args.port}/"
    try:
        server, _ = make_server(args.port)
    except OSError as e:
        if _already_running(url):  # started twice: show the window that is there
            if not args.no_browser:
                webbrowser.open(url)
            return 0
        print(f"Port {args.port} cannot be used ({e}). Try another one with --port.")
        return 1
    print(f"Foci Analyzer {__version__} at {url}  (Ctrl+C to stop)")
    if not args.no_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        models.stop_all()
    return 0
