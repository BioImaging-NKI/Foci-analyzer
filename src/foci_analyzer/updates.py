"""Updates: is there a newer version on GitHub, and installing it.

The check asks the GitHub API for the releases of the repository (nothing else is sent) at most once a day,
in the background, and on the Settings page. The repository also holds the releases of the Fiji macro, so
only releases with the Windows setup attached (`FociAnalyzer-Setup-*.exe`) count.

The update downloads the release's source archive into the installation folder and runs its
`installer/install.py --update` detached, with the networks that are installed now. That script waits until
this program has stopped, shows its progress on the window's port (the browser tab follows it), installs,
and starts the new version. The user's settings, experiment and results are not touched.

Only an installation made by the setup or install.sh can update itself (it has `bin/uv` and the program
folder); a developer's checkout is updated with git.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import urllib.request
import zipfile
from pathlib import Path

from . import __version__

REPOSITORY = "BioImaging-NKI/Foci-analyzer"
ASSET = re.compile(r"FociAnalyzer-Setup-.+\.exe$")
DAY = 24 * 3600


def releases_url() -> str:
    # FOCI_ANALYZER_UPDATE_URL: for tests, a file:// or http:// URL with the same JSON as the GitHub API
    return os.environ.get("FOCI_ANALYZER_UPDATE_URL") or f"https://api.github.com/repos/{REPOSITORY}/releases?per_page=30"


def version_key(text: str) -> tuple:
    """'v1.2.3' -> (1, 2, 3, 1, 0); a development version sorts before its release: '1.2.3.dev4' -> (1, 2, 3, 0, 4)."""
    m = re.match(r"v?(\d+)\.(\d+)(?:\.(\d+))?(?:\.?(dev|a|b|rc)(\d*))?", str(text).strip())
    if not m:
        return (0,)
    return (int(m[1]), int(m[2]), int(m[3] or 0), 0 if m[4] else 1, int(m[5] or 0))


def latest(timeout: float = 15) -> dict | None:
    """The newest Foci Analyzer release (not a draft, not a pre-release, with the setup attached), or None."""
    req = urllib.request.Request(releases_url(), headers={"Accept": "application/vnd.github+json",
                                                          "User-Agent": "Foci-Analyzer-update-check"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        found = json.loads(r.read().decode("utf-8"))
    best = None
    for rel in found if isinstance(found, list) else []:
        if rel.get("draft") or rel.get("prerelease"):
            continue
        setup = next((a.get("browser_download_url") for a in rel.get("assets", []) if ASSET.search(a.get("name", ""))), None)
        if not setup:
            continue  # a release of the Fiji macro
        info = {"version": str(rel.get("tag_name", "")).lstrip("v"), "notes": rel.get("body") or "",
                "page": rel.get("html_url", ""), "setup": setup, "source": rel.get("zipball_url", ""),
                "date": (rel.get("published_at") or "")[:10]}
        if best is None or version_key(info["version"]) > version_key(best["version"]):
            best = info
    return best


def is_newer(version: str) -> bool:
    return version_key(version) > version_key(__version__)


# ---- installing

def layout() -> dict | None:
    """The installation this program runs from, when it can update itself: its folder, uv and the
    program's launcher. None for a developer's checkout."""
    target = Path(sys.prefix).resolve().parent
    win = sys.platform == "win32"
    uv = target / "bin" / ("uv.exe" if win else "uv")
    if not uv.is_file() or not (target / "python").is_dir():
        return None
    return {"target": target, "uv": uv, "windows": win}


def networks_installed(target: Path, installed: dict, python_for) -> list[str]:
    """The networks installed in this folder (not those of the user's own environments): installed again."""
    names = []
    for op in ("stardist", "cellpose", "cellpose4", "instanseg", "spotiflow", "trackastra"):
        py = python_for(op)
        if installed.get(op) and py and Path(py).resolve().is_relative_to(target):
            names.append(op)
    if "cellpose" in names and not installed.get("cellpose3"):
        names.remove("cellpose")  # that Cellpose is a Cellpose 4 (cellpose4-env): installed as cellpose4
    return names


def download(url: str, dest: Path, progress) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "Foci-Analyzer-update"})
    with urllib.request.urlopen(req, timeout=60) as r, dest.open("wb") as f:
        total = int(r.headers.get("Content-Length") or 0)
        done = 0
        while chunk := r.read(1 << 16):
            f.write(chunk)
            done += len(chunk)
            progress(done, total)


def start(info: dict, port: int, networks: list[str], progress) -> None:
    """Download and unpack the new version, then hand over to its install.py (detached). The caller stops
    this program right after."""
    where = layout()
    if where is None:
        raise RuntimeError("this Foci Analyzer was not installed with the setup or install.sh; update it with git")
    work = where["target"] / "update"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    archive = work / "source.zip"
    download(info["source"], archive, progress)
    with zipfile.ZipFile(archive) as z:
        z.extractall(work / "source")
    script = next((work / "source").glob("*/installer/install.py"), None) or work / "source" / "installer" / "install.py"
    if not script.is_file():
        raise RuntimeError("the downloaded version has no installer/install.py")
    target = where["target"]
    cmd = [str(where["uv"]), "run", "--no-project", "--python", "3.12", str(script), "--target", str(target),
           "--uv", str(where["uv"]), "--networks", ",".join(networks) or "none", "--update",
           "--wait", str(os.getpid()), "--serve", str(port)]
    if not where["windows"]:
        cmd.append("--shortcut")
    env = {**os.environ, "UV_PYTHON_INSTALL_DIR": str(target / "python"), "UV_CACHE_DIR": str(target / "cache"),
           "UV_NO_PROGRESS": "1", "PYTHONIOENCODING": "utf-8"}
    env.pop("VIRTUAL_ENV", None)
    out = (work / "update-output.log").open("wb")
    if where["windows"]:
        flags = subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW
        subprocess.Popen(cmd, cwd=str(work), env=env, stdout=out, stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, creationflags=flags, close_fds=True)
    else:
        subprocess.Popen(cmd, cwd=str(work), env=env, stdout=out, stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, start_new_session=True, close_fds=True)


# ---- what the user chose (in the Foci Analyzer folder, next to program.json)

class Preferences:
    def __init__(self, path: Path):
        self.path = path
        self.data = {"check": True, "last_check": 0, "snooze_until": 0, "skip": "", "latest": None}
        try:
            self.data.update(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, ValueError):
            pass
        self.lock = threading.Lock()

    def save(self) -> None:
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.write_text(json.dumps(self.data, indent=1), encoding="utf-8")
        except OSError:
            pass

    def due(self) -> bool:
        return bool(self.data["check"]) and time.time() - float(self.data.get("last_check") or 0) > DAY

    def show(self) -> bool:
        """Whether the window should offer the newer version now (not skipped, not snoozed)."""
        new = self.data.get("latest")
        return bool(new and is_newer(new["version"]) and new["version"] != self.data.get("skip")
                    and time.time() >= float(self.data.get("snooze_until") or 0))
