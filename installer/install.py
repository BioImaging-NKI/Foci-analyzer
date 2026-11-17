"""Install Foci Analyzer and the neural networks; the same steps on Windows, macOS and Linux.

Run by the installers (setup.exe on Windows, install.sh on macOS and Linux, install.bat for developers) as
    uv run --no-project --python 3.12 installer/install.py [options]
It needs only the standard library and uv. Without --networks and with a terminal it asks its questions.

Everything goes into one folder (--target): the program environment (.venv), the network environments
(models-env: PyTorch networks; cellpose4-env: Cellpose 4, which cannot live next to Cellpose 3; stardist-env:
StarDist, which needs TensorFlow), the Python versions and the
download cache. Removing that folder removes the program. The program finds the network environments next
to its own; environments of the user are written to program.json in the user's Foci Analyzer folder.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent  # the program: pyproject.toml, uv.lock, src
NETWORKS = ("stardist", "cellpose", "cellpose4", "instanseg", "spotiflow", "trackastra")
TORCH_NETWORKS = ("cellpose", "instanseg", "spotiflow", "trackastra")
WINDOWS, MAC = sys.platform == "win32", sys.platform == "darwin"
# versions tested together (plan.md, section 9); a newer TensorFlow is a hundred times slower on Windows
TORCH = ["torch==2.6.0", "torchvision==0.21.0"]
PACKAGES = {"cellpose": ["cellpose==3.1.1.3"], "instanseg": ["instanseg-torch==0.1.2"],
            "spotiflow": ["spotiflow==0.6.5", "csbdeep"], "trackastra": ["trackastra==0.5.6"]}
CELLPOSE4 = ["cellpose==4.2.1.1"]  # its own environment: Cellpose 3 and 4 cannot be installed together
STARDIST = ["stardist==0.9.2", "csbdeep==0.8.2", "numpy<2", "scikit-image"]
IMPORTS = {"stardist": "stardist", "cellpose": "cellpose", "cellpose4": "cellpose", "instanseg": "instanseg", "spotiflow": "spotiflow",
           "trackastra": "trackastra"}
LOG: Path | None = None
# disk space in GB once installed (measured on Windows; the setup shows the same numbers)
SIZE_PROGRAM = 0.5
SIZE_TORCH = {"gpu": 4.4, "cpu": 1.0}
SIZE = {"cellpose": 0.1, "instanseg": 0.3, "spotiflow": 0.2, "trackastra": 0.2, "stardist": 1.6, "cellpose4": 0.1}


class Progress:
    """For the Windows setup (--progress FILE): how far the installation is, the current step and the messages,
    in a small text file that the setup window reads every half second. Line 1: per mille done; line 2: the
    step; line 3: the number of messages so far; then the last messages.

    Steps count by their size. Within a step uv tells what it downloads ("Downloading torch (2.4GiB)", then
    "Downloaded torch" when done, as run() writes it); a download in progress is estimated from the time it takes."""

    SPEED = 8e6  # bytes per second assumed for a download in progress, until one finished

    def __init__(self, path: Path | None):
        self.path, self.lines, self.count = path, [], 0
        self.weights, self.step, self.status = [], -1, "Starting"
        self.files: dict[str, float] = {}  # file being downloaded -> bytes
        self.done_bytes = self.total_bytes = 0.0
        self.since, self.speed = time.time(), self.SPEED
        self.final = False
        if path:
            import threading
            threading.Thread(target=self._tick, daemon=True).start()

    def plan(self, weights: list[float]) -> None:
        self.weights = [max(w, 0.01) for w in weights]

    def next(self, status: str) -> None:
        self.step += 1
        self.status = f"Step {self.step + 1} of {len(self.weights)}: {status}"
        self.files.clear()
        self.done_bytes = self.total_bytes = 0.0
        self.since = time.time()
        self.write()

    def line(self, text: str) -> None:
        self.lines.append(text)
        self.count += 1
        del self.lines[:-300]
        m = re.match(r"Downloading (\S+) \(([\d.]+)\s*([KMG])iB\)", text)
        if m:
            size = float(m[2]) * {"K": 2**10, "M": 2**20, "G": 2**30}[m[3]]
            self.files[m[1]] = size
            self.total_bytes += size
            self.since = time.time()
        m = re.match(r"\s+Downloaded (\S+)$", text)
        if m and m[1] in self.files:
            size = self.files.pop(m[1])
            self.done_bytes += size
            seconds = time.time() - self.since
            if size > 50e6 and seconds > 2:
                self.speed = max(self.speed, size / seconds)

    def permille(self) -> int:
        if self.final:
            return 1000
        if not self.weights or self.step < 0:
            return 0
        part = 0.0
        if self.total_bytes:
            busy = sum(self.files.values())
            guess = min(0.95 * busy, (time.time() - self.since) * self.speed)
            part = 0.85 * (self.done_bytes + guess) / self.total_bytes
        if not self.files and self.total_bytes and self.done_bytes >= self.total_bytes:
            part = 0.9  # downloaded; installing
        total = sum(self.weights)
        return int(1000 * (sum(self.weights[:self.step]) + part * self.weights[self.step]) / total)

    def write(self) -> None:
        if not self.path:
            return
        text = "\n".join([str(self.permille()), self.status, str(self.count), *self.lines]) + "\n"
        tmp = self.path.with_suffix(".tmp")
        try:
            tmp.write_text(text, encoding="ascii", errors="replace")
            os.replace(tmp, self.path)
        except OSError:
            pass  # the setup is reading it: next time

    def done(self) -> None:
        self.final, self.status = True, "Ready"
        self.write()

    def _tick(self) -> None:
        while True:
            time.sleep(0.5)
            self.write()


PROGRESS = Progress(None)


def say(text: str) -> None:
    print(text, flush=True)
    for line in text.strip("\n").splitlines() or [""]:
        PROGRESS.line(line)
    if LOG:
        with LOG.open("a", encoding="utf-8") as f:
            f.write(f"{time.strftime('%H:%M:%S')} {text}\n")


def run(uv: str, *args: str, env: dict | None = None) -> None:
    cmd = [uv, *args]
    if LOG:
        with LOG.open("a", encoding="utf-8") as f:
            f.write("$ " + " ".join(cmd) + "\n")
    # uv's messages go through say(): to the console, the log and the setup window
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env={**os.environ, **(env or {})},
                         text=True, encoding="utf-8", errors="replace", bufsize=1)
    for line in p.stdout:
        line = re.sub(r"\x1b\[[0-9;]*m", "", line.rstrip())
        say(re.sub(r"^(\s+)Downloading (\S+)$", r"\1Downloaded \2", line))  # uv's "done" line reads oddly
    if p.wait():
        raise SystemExit(f"\nThis step failed (exit code {p.returncode}): {' '.join(cmd)}\n"
                         "The messages above say why; an internet connection is needed.")


def python_of(env: Path) -> Path:
    return env / ("Scripts/python.exe" if WINDOWS else "bin/python")


def venv(uv: str, path: Path, python: str, env: dict) -> None:
    """An environment, unless it is there already (its python.exe may be in use by a running program). One that
    is only partly there (an uninstall while the program ran) is removed and made again."""
    if path.exists() and not healthy(path):
        say(f"Removing the incomplete environment {path}")
        shutil.rmtree(path, ignore_errors=True)
    if not python_of(path).is_file():
        run(uv, "venv", str(path), "--python", python, "--no-project", env=env)


def healthy(path: Path) -> bool:
    """Whether this environment's Python starts."""
    if not (python_of(path).is_file() and (path / "pyvenv.cfg").is_file()):
        return False
    try:
        return subprocess.run([str(python_of(path)), "-c", "import os, json"], check=False, capture_output=True,
                              timeout=60).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def nvidia() -> bool:
    if MAC:
        return False
    try:
        return subprocess.run(["nvidia-smi", "-L"], check=False, capture_output=True, timeout=20).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def conda_env(python: str) -> dict:
    """Environment variables for a Python of a conda environment that is not activated (the same as
    models.env_for in the program, which this script cannot import)."""
    env = dict(os.environ)
    root = Path(python).parent
    if not (root / "conda-meta").is_dir():
        root = root.parent
    if (root / "conda-meta").is_dir():
        dirs = ([root, root / "Library" / "mingw-w64" / "bin", root / "Library" / "usr" / "bin", root / "Library" / "bin",
                 root / "Scripts", root / "bin"] if WINDOWS else [root / "bin"])
        env["PATH"] = os.pathsep.join([*map(str, dirs), env.get("PATH", "")])
        env["CONDA_PREFIX"] = str(root)
    return env


def check_python(path: str, network: str) -> str:
    """The version of the network in another environment, or why it cannot be used."""
    if not Path(path).is_file():
        return f"not found: {path}"
    mod = IMPORTS[network]
    dist = "instanseg-torch" if network == "instanseg" else mod
    code = f"import importlib.metadata as m, importlib; importlib.import_module('{mod}'); print(m.version('{dist}'))"
    try:
        r = subprocess.run([path, "-c", code], check=False, capture_output=True, text=True, timeout=300,
                           env=conda_env(path))
    except (OSError, subprocess.TimeoutExpired) as e:
        return f"cannot run it: {e}"
    last = r.stdout.strip().splitlines()[-1:] or ["?"]  # the networks print their own messages first
    return f"version {last[0]}" if r.returncode == 0 else f"{network} is not in this environment"


def ask(args) -> None:
    """The questions, in a terminal, when the installer did not get them as options."""
    print("\nNeural networks (each can be left out; the program works without them, more slowly or with a")
    print("simpler method):")
    print("  StarDist (nuclei; recommended), Cellpose 3 (nuclei or cells, 2D and 3D), Cellpose 4 (the newest")
    print("  Cellpose models, slower), InstanSeg (nuclei or cells), Spotiflow (foci), Trackastra (tracking of cells")
    print("  in time-lapses). Together about 10 GB.")
    a = input("Install all of them? [Y/n/choose] ").strip().lower()
    if a in ("", "y", "yes"):
        args.networks = ",".join(NETWORKS)
    elif a.startswith("c"):
        args.networks = ",".join(n for n in NETWORKS if input(f"  {n}? [Y/n] ").strip().lower() in ("", "y", "yes"))
    else:
        args.networks = "none"
    if input("Use networks from an environment you already have (for example your own Cellpose)? [y/N] ").strip().lower() in ("y", "yes"):
        for n in NETWORKS:
            p = input(f"  python of the environment with {n} (empty = none): ").strip().strip('"')
            if p:
                args.own.append(f"{n}={p}")


def _pe_end(data: bytes) -> int:
    """Where the program part of an .exe ends (what follows is appended data, such as uv's start script)."""
    import struct
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    sections, opt_size = struct.unpack_from("<H", data, pe + 6)[0], struct.unpack_from("<H", data, pe + 20)[0]
    table = pe + 24 + opt_size
    end = max(struct.unpack_from("<II", data, table + 40 * i + 16)[1] + struct.unpack_from("<II", data, table + 40 * i + 16)[0]
              for i in range(sections))
    return end


def branded_launcher(src: Path, dst: Path, ico: Path) -> bool:
    """A copy of uv's launcher `src` with the icon of `ico` (Windows). uv's launcher finds its Python and script
    in data appended after the program part, so that part is set aside while the icon is replaced, and appended
    again. Returns False (and leaves a plain copy) when the icon cannot be set."""
    import ctypes
    import struct
    data = src.read_bytes()
    end = _pe_end(data)
    dst.write_bytes(data)
    if not ico.is_file() or end >= len(data):
        return False
    try:
        icon = ico.read_bytes()
        _, kind, count = struct.unpack_from("<HHH", icon, 0)
        if kind != 1:
            return False
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.BeginUpdateResourceW.restype = ctypes.c_void_p
        k32.BeginUpdateResourceW.argtypes = [ctypes.c_wchar_p, ctypes.c_bool]
        k32.UpdateResourceW.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_ushort,
                                        ctypes.c_char_p, ctypes.c_uint]
        k32.EndUpdateResourceW.argtypes = [ctypes.c_void_p, ctypes.c_bool]
        dst.write_bytes(data[:end])  # the program part only; the rest is appended again below
        h = k32.BeginUpdateResourceW(str(dst), False)
        if not h:
            raise OSError(ctypes.get_last_error())
        group = struct.pack("<HHH", 0, 1, count)
        for i in range(count):
            w, hh, colours, res, planes, bits, size, offset = struct.unpack_from("<BBBBHHII", icon, 6 + 16 * i)
            k32.UpdateResourceW(h, 3, 100 + i, 0, icon[offset:offset + size], size)  # RT_ICON
            group += struct.pack("<BBBBHHIH", w, hh, colours, res, planes, bits, size, 100 + i)
        for name in (1, 101):  # the group icon that Windows shows: the first one, by number
            k32.UpdateResourceW(h, 14, name, 0, group, len(group))  # RT_GROUP_ICON
        if not k32.EndUpdateResourceW(h, False):
            raise OSError(ctypes.get_last_error())
        with dst.open("ab") as f:
            f.write(data[end:])
        return True
    except (OSError, ValueError, AttributeError, struct.error) as e:
        say(f"  (the program icon could not be set: {e})")
        dst.write_bytes(data)
        return False


def shortcuts(target: Path, desktop: bool, menu: bool = True) -> None:
    exe = python_of(target / ".venv").parent / ("foci-analyzer-window.exe" if WINDOWS else "foci-analyzer-window")
    icons = target / "icons"  # kept with the program: the folder the installer ran from may be removed later
    icons.mkdir(exist_ok=True)
    for name in ("FociAnalyzer.ico", "FociAnalyzer.icns", "logo-512.png"):
        if (HERE / name).is_file():
            shutil.copy2(HERE / name, icons / name)
    if WINDOWS and exe.is_file():  # "Foci Analyzer.exe" in the program folder, with the logo as its icon
        launcher = target / "Foci Analyzer.exe"
        if branded_launcher(exe, launcher, icons / "FociAnalyzer.ico"):
            say(f"Made {launcher}")
        exe = launcher
    if not menu:
        return
    if WINDOWS:
        places = ["[Environment]::GetFolderPath('Programs')"] + (["[Environment]::GetFolderPath('Desktop')"] if desktop else [])
        for place in places:
            ps = (f"$s = (New-Object -ComObject WScript.Shell).CreateShortcut([IO.Path]::Combine({place}, 'Foci Analyzer.lnk')); "
                  f"$s.TargetPath = '{exe}'; $s.WorkingDirectory = '{target}'; $s.Description = 'Foci Analyzer'; "
                  f"$s.IconLocation = '{icons / 'FociAnalyzer.ico'}'; $s.Save()")
            subprocess.run(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-Command", ps], check=False)
    elif MAC:  # a small app in ~/Applications that starts the program (made here, so macOS does not block it)
        app = Path.home() / "Applications" / "Foci Analyzer.app" / "Contents"
        (app / "MacOS").mkdir(parents=True, exist_ok=True)
        (app / "Resources").mkdir(exist_ok=True)
        if (icons / "FociAnalyzer.icns").is_file():
            shutil.copy2(icons / "FociAnalyzer.icns", app / "Resources" / "FociAnalyzer.icns")
        start = app / "MacOS" / "Foci Analyzer"
        start.write_text(f'#!/bin/sh\nexec "{exe}"\n', encoding="utf-8")
        start.chmod(0o755)
        (app / "Info.plist").write_text(
            '<?xml version="1.0" encoding="UTF-8"?>\n<plist version="1.0"><dict>'
            "<key>CFBundleName</key><string>Foci Analyzer</string>"
            "<key>CFBundleExecutable</key><string>Foci Analyzer</string>"
            "<key>CFBundleIconFile</key><string>FociAnalyzer</string>"
            "<key>CFBundleIdentifier</key><string>nl.nki.foci-analyzer</string>"
            "<key>CFBundlePackageType</key><string>APPL</string></dict></plist>\n", encoding="utf-8")
        say(f"Made {app.parent}")
    else:
        share = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))
        apps = share / "applications"
        apps.mkdir(parents=True, exist_ok=True)
        icon = str(icons / "logo-512.png")
        if (icons / "logo-512.png").is_file():  # in the icon theme too, where every desktop looks for it
            themed = share / "icons" / "hicolor" / "512x512" / "apps"
            themed.mkdir(parents=True, exist_ok=True)
            shutil.copy2(icons / "logo-512.png", themed / "foci-analyzer.png")
            icon = "foci-analyzer"
        (apps / "foci-analyzer.desktop").write_text(
            "[Desktop Entry]\nType=Application\nName=Foci Analyzer\nComment=Foci in nuclei and cells\n"
            f'Exec="{exe}"\nIcon={icon}\nTerminal=false\nCategories=Science;Education;\n'
            "StartupWMClass=Foci Analyzer\n", encoding="utf-8")
        say(f"Made {apps / 'foci-analyzer.desktop'}")


def main(argv: list[str] | None = None) -> None:
    global LOG, PROGRESS
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--target", required=True, help="the folder to install in")
    p.add_argument("--uv", default=os.environ.get("UV") or shutil.which("uv") or "uv", help="the uv to use")
    p.add_argument("--networks", help=f"comma-separated: {','.join(NETWORKS)}; 'all' or 'none'")
    p.add_argument("--own", action="append", default=[], metavar="NETWORK=PYTHON",
                   help="use this network from another environment (its python); may be repeated")
    p.add_argument("--device", choices=["auto", "gpu", "cpu"], default="auto", help="PyTorch for the NVIDIA card or the processor")
    p.add_argument("--home", default=os.environ.get("FOCI_ANALYZER_HOME") or str(Path.home() / ".foci_analyzer"),
                   help="the user's Foci Analyzer folder (for program.json)")
    p.add_argument("--networks-only", action="store_true", help="only install the networks (the program is there already)")
    p.add_argument("--shortcut", action="store_true", help="add the program to the start menu / Applications")
    p.add_argument("--desktop", action="store_true", help="also put a shortcut on the desktop (Windows)")
    p.add_argument("--progress", help="a file in which the Windows setup follows the installation")
    args = p.parse_args(argv)
    if args.networks is None:
        if sys.stdin.isatty():
            ask(args)
        else:
            args.networks = "none"

    target = Path(args.target).resolve()
    target.mkdir(parents=True, exist_ok=True)
    LOG = target / "install.log"
    PROGRESS = Progress(Path(args.progress) if args.progress else None)
    own = {}
    for item in args.own:
        name, _, path = item.partition("=")
        if name not in NETWORKS or not path.strip():
            raise SystemExit(f"--own needs NETWORK=PYTHON with NETWORK one of {', '.join(NETWORKS)}: {item}")
        own[name] = path.strip().strip('"')
    wanted = set(NETWORKS) if args.networks == "all" else set() if args.networks in ("none", "") else set(args.networks.split(","))
    if wanted - set(NETWORKS):
        raise SystemExit(f"unknown network(s): {', '.join(sorted(wanted - set(NETWORKS)))}")
    wanted -= set(own)  # a network from the user's own environment is not installed again
    # Python versions and the download cache inside the program folder: removing the folder removes all.
    # (The cache costs little extra disk space: uv links the files from it into the environments.)
    env = {"UV_PYTHON_INSTALL_DIR": str(target / "python"), "UV_CACHE_DIR": str(target / "cache"),
           "UV_PYTHON_PREFERENCE": "only-managed"}
    gpu = nvidia() if args.device == "auto" else args.device == "gpu"
    torch_size = SIZE_TORCH["gpu" if gpu else "cpu"]
    torch_nets = sorted(wanted & set(TORCH_NETWORKS))
    weights = ([] if args.networks_only else [SIZE_PROGRAM]) + (
        [torch_size + sum(SIZE[n] for n in torch_nets)] if torch_nets else []) + (
        [torch_size + SIZE["cellpose4"]] if "cellpose4" in wanted else []) + (
        [SIZE["stardist"]] if "stardist" in wanted else []) + ([0.01] if own else []) + (
        [] if args.networks_only else [0.01])
    PROGRESS.plan(weights)
    steps, step = len(weights), 0

    def next_step(text, short):
        nonlocal step
        step += 1
        PROGRESS.next(short)
        say(f"\n[{step}/{steps}] {text}")

    say(f"Foci Analyzer: installing in {target} ({platform.system()} {platform.machine()})")
    if not args.networks_only:
        next_step("The program", "Python and the program")
        if (target / ".venv").exists() and not healthy(target / ".venv"):
            say("Removing the incomplete program environment")
            shutil.rmtree(target / ".venv", ignore_errors=True)
        run(args.uv, "sync", "--project", str(SOURCE), "--no-dev", "--frozen", "--no-editable", "--python", "3.12",
            env={**env, "UV_PROJECT_ENVIRONMENT": str(target / ".venv")})
    index = [] if MAC else ["--index-url", "https://download.pytorch.org/whl/" + ("cu124" if gpu else "cpu"),
                            "--extra-index-url", "https://pypi.org/simple", "--index-strategy", "unsafe-best-match"]

    if torch_nets:
        where = ("the Apple graphics chip (MPS) where it can" if MAC else
                 "the NVIDIA graphics card" if gpu else "the processor (no NVIDIA graphics card found)")
        next_step(f"PyTorch networks ({', '.join(torch_nets)}), running on {where}; this takes a while",
                  "PyTorch and " + ", ".join(n.capitalize() for n in torch_nets))
        menv = target / "models-env"
        venv(args.uv, menv, "3.12", env)
        run(args.uv, "pip", "install", "--python", str(python_of(menv)), *TORCH, *index, env=env)
        run(args.uv, "pip", "install", "--python", str(python_of(menv)), *TORCH,
            *[pkg for n in torch_nets for pkg in PACKAGES[n]], *index, env=env)

    if "cellpose4" in wanted:
        next_step("Cellpose 4 (its own environment); this takes a while", "Cellpose 4 (with its own PyTorch)")
        cenv = target / "cellpose4-env"
        venv(args.uv, cenv, "3.12", env)
        run(args.uv, "pip", "install", "--python", str(python_of(cenv)), *TORCH, *CELLPOSE4, *index, env=env)

    if "stardist" in wanted:
        gpu_tf = not WINDOWS and not MAC and nvidia() and args.device != "cpu"
        next_step("StarDist (TensorFlow" + (", NVIDIA graphics card)" if gpu_tf else ", processor)"), "StarDist (TensorFlow)")
        senv = target / "stardist-env"
        venv(args.uv, senv, "3.11", env)
        tf = "tensorflow[and-cuda]==2.15.1" if gpu_tf else "tensorflow==2.15.1"
        run(args.uv, "pip", "install", "--python", str(python_of(senv)), tf, *STARDIST, env=env)

    if own:
        next_step("Your own environments", "Your own environments")
        home = Path(args.home)
        home.mkdir(parents=True, exist_ok=True)
        pj = home / "program.json"
        try:
            program = json.loads(pj.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            program = {}
        program.setdefault("python", "")
        program.setdefault("pythons", {})
        for name, path in own.items():
            say(f"  {name}: {path} -> {check_python(path, name)}")
            program["pythons"][name] = path
        pj.write_text(json.dumps(program, indent=1), encoding="utf-8")

    if args.networks_only:
        say("\nReady.")
        PROGRESS.done()
        return
    next_step("Shortcuts", "Shortcuts")
    shortcuts(target, args.desktop, menu=args.shortcut)  # on Windows always "Foci Analyzer.exe" with the logo
    exe = target / "Foci Analyzer.exe" if WINDOWS else python_of(target / ".venv").parent / "foci-analyzer"
    say(f"\nReady. Start Foci Analyzer {'from the start menu' if args.shortcut and WINDOWS else 'from Applications' if args.shortcut and MAC else 'from the applications menu' if args.shortcut else 'with'}"
        f"{'' if args.shortcut else ' ' + str(exe)}; it opens in the browser.")
    PROGRESS.done()


if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        if e.code not in (0, None):
            say(str(e.code))
            PROGRESS.status = "The installation did not finish: " + str(e.code).strip().splitlines()[0]
            PROGRESS.write()
        raise
