"""Install Foci Analyzer and the neural networks; the same steps on Windows, macOS and Linux.

Run by the installers (setup.exe on Windows, install.sh on macOS and Linux, install.bat for developers) as
    uv run --no-project --python 3.12 installer/install.py [options]
It needs only the standard library and uv. Without --networks and with a terminal it asks its questions.

Everything goes into one folder (--target): the program environment (.venv), the network environments
(models-env: PyTorch networks; stardist-env: StarDist, which needs TensorFlow), the Python versions and the
download cache. Removing that folder removes the program. The program finds the network environments next
to its own; environments of the user are written to program.json in the user's Foci Analyzer folder.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent  # the program: pyproject.toml, uv.lock, src
NETWORKS = ("stardist", "cellpose", "instanseg", "spotiflow", "trackastra")
TORCH_NETWORKS = ("cellpose", "instanseg", "spotiflow", "trackastra")
WINDOWS, MAC = sys.platform == "win32", sys.platform == "darwin"
# versions tested together (plan.md, section 9); a newer TensorFlow is a hundred times slower on Windows
TORCH = ["torch==2.6.0", "torchvision==0.21.0"]
PACKAGES = {"cellpose": ["cellpose==3.1.1.3"], "instanseg": ["instanseg-torch==0.1.2"],
            "spotiflow": ["spotiflow==0.6.5", "csbdeep"], "trackastra": ["trackastra==0.5.6"]}
STARDIST = ["stardist==0.9.2", "csbdeep==0.8.2", "numpy<2", "scikit-image"]
IMPORTS = {"stardist": "stardist", "cellpose": "cellpose", "instanseg": "instanseg", "spotiflow": "spotiflow",
           "trackastra": "trackastra"}
LOG: Path | None = None


def say(text: str) -> None:
    print(text, flush=True)
    if LOG:
        with LOG.open("a", encoding="utf-8") as f:
            f.write(f"{time.strftime('%H:%M:%S')} {text}\n")


def run(uv: str, *args: str, env: dict | None = None) -> None:
    cmd = [uv, *args]
    if LOG:
        with LOG.open("a", encoding="utf-8") as f:
            f.write("$ " + " ".join(cmd) + "\n")
    r = subprocess.run(cmd, check=False, env={**os.environ, **(env or {})})
    if r.returncode:
        raise SystemExit(f"\nThis step failed (exit code {r.returncode}): {' '.join(cmd)}\n"
                         "The messages above say why; an internet connection is needed.")


def python_of(env: Path) -> Path:
    return env / ("Scripts/python.exe" if WINDOWS else "bin/python")


def nvidia() -> bool:
    if MAC:
        return False
    try:
        return subprocess.run(["nvidia-smi", "-L"], check=False, capture_output=True, timeout=20).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def check_python(path: str, network: str) -> str:
    """The version of the network in another environment, or why it cannot be used."""
    if not Path(path).is_file():
        return f"not found: {path}"
    mod = IMPORTS[network]
    dist = "instanseg-torch" if network == "instanseg" else network
    code = f"import importlib.metadata as m, importlib; importlib.import_module('{mod}'); print(m.version('{dist}'))"
    try:
        r = subprocess.run([path, "-c", code], check=False, capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as e:
        return f"cannot run it: {e}"
    last = r.stdout.strip().splitlines()[-1:] or ["?"]  # the networks print their own messages first
    return f"version {last[0]}" if r.returncode == 0 else f"{network} is not in this environment"


def ask(args) -> None:
    """The questions, in a terminal, when the installer did not get them as options."""
    print("\nNeural networks (each can be left out; the program works without them, more slowly or with a")
    print("simpler method):")
    print("  StarDist (nuclei; recommended), Cellpose (nuclei or cells, 2D and 3D), InstanSeg (nuclei or cells),")
    print("  Spotiflow (foci), Trackastra (tracking of cells in time-lapses). Together about 7 GB.")
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


def shortcuts(target: Path, desktop: bool) -> None:
    exe = python_of(target / ".venv").parent / ("foci-analyzer-window.exe" if WINDOWS else "foci-analyzer-window")
    icons = target / "icons"  # kept with the program: the folder the installer ran from may be removed later
    icons.mkdir(exist_ok=True)
    for name in ("FociAnalyzer.ico", "FociAnalyzer.icns", "logo-512.png"):
        if (HERE / name).is_file():
            shutil.copy2(HERE / name, icons / name)
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
        apps = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share")) / "applications"
        apps.mkdir(parents=True, exist_ok=True)
        (apps / "foci-analyzer.desktop").write_text(
            "[Desktop Entry]\nType=Application\nName=Foci Analyzer\nComment=Foci in nuclei and cells\n"
            f'Exec="{exe}"\nIcon={icons / "logo-512.png"}\nTerminal=false\nCategories=Science;Education;\n', encoding="utf-8")
        say(f"Made {apps / 'foci-analyzer.desktop'}")


def main(argv: list[str] | None = None) -> None:
    global LOG
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--target", required=True, help="the folder to install in")
    p.add_argument("--uv", default=os.environ.get("UV") or shutil.which("uv") or "uv", help="the uv to use")
    p.add_argument("--networks", help=f"comma-separated: {','.join(NETWORKS)}; 'all' or 'none'")
    p.add_argument("--own", action="append", default=[], metavar="NETWORK=PYTHON",
                   help="use this network from another environment (its python); may be repeated")
    p.add_argument("--device", choices=["auto", "gpu", "cpu"], default="auto", help="PyTorch for the NVIDIA card or the processor")
    p.add_argument("--home", default=os.environ.get("FOCI_ANALYZER_HOME") or str(Path.home() / ".foci_analyzer"),
                   help="the user's Foci Analyzer folder (for program.json)")
    p.add_argument("--shortcut", action="store_true", help="add the program to the start menu / Applications")
    p.add_argument("--desktop", action="store_true", help="also put a shortcut on the desktop (Windows)")
    args = p.parse_args(argv)
    if args.networks is None:
        if sys.stdin.isatty():
            ask(args)
        else:
            args.networks = "none"

    target = Path(args.target).resolve()
    target.mkdir(parents=True, exist_ok=True)
    LOG = target / "install.log"
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
    steps = 1 + bool(wanted & set(TORCH_NETWORKS)) + ("stardist" in wanted) + bool(own) + 1
    step = 0

    def next_step(text):
        nonlocal step
        step += 1
        say(f"\n[{step}/{steps}] {text}")

    say(f"Foci Analyzer: installing in {target} ({platform.system()} {platform.machine()})")
    next_step("The program")
    run(args.uv, "sync", "--project", str(SOURCE), "--no-dev", "--frozen", "--no-editable", "--python", "3.12",
        env={**env, "UV_PROJECT_ENVIRONMENT": str(target / ".venv")})

    torch_nets = sorted(wanted & set(TORCH_NETWORKS))
    if torch_nets:
        gpu = nvidia() if args.device == "auto" else args.device == "gpu"
        where = ("the Apple graphics chip (MPS) where it can" if MAC else
                 "the NVIDIA graphics card" if gpu else "the processor (no NVIDIA graphics card found)")
        next_step(f"PyTorch networks ({', '.join(torch_nets)}), running on {where}; this takes a while")
        menv = target / "models-env"
        run(args.uv, "venv", str(menv), "--python", "3.12", "--allow-existing", "--no-project", env=env)
        index = [] if MAC else ["--index-url", "https://download.pytorch.org/whl/" + ("cu124" if gpu else "cpu"),
                                "--extra-index-url", "https://pypi.org/simple", "--index-strategy", "unsafe-best-match"]
        run(args.uv, "pip", "install", "--python", str(python_of(menv)), *TORCH, *index, env=env)
        run(args.uv, "pip", "install", "--python", str(python_of(menv)), *TORCH,
            *[pkg for n in torch_nets for pkg in PACKAGES[n]], *index, env=env)

    if "stardist" in wanted:
        gpu_tf = not WINDOWS and not MAC and nvidia() and args.device != "cpu"
        next_step("StarDist (TensorFlow" + (", NVIDIA graphics card)" if gpu_tf else ", processor)"))
        senv = target / "stardist-env"
        run(args.uv, "venv", str(senv), "--python", "3.11", "--allow-existing", "--no-project", env=env)
        tf = "tensorflow[and-cuda]==2.15.1" if gpu_tf else "tensorflow==2.15.1"
        run(args.uv, "pip", "install", "--python", str(python_of(senv)), tf, *STARDIST, env=env)

    if own:
        next_step("Your own environments")
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

    next_step("Shortcuts")
    if args.shortcut:
        shortcuts(target, args.desktop)
    exe = python_of(target / ".venv").parent / ("foci-analyzer-window.exe" if WINDOWS else "foci-analyzer")
    say(f"\nReady. Start Foci Analyzer {'from the start menu' if args.shortcut and WINDOWS else 'from Applications' if args.shortcut and MAC else 'from the applications menu' if args.shortcut else 'with'}"
        f"{'' if args.shortcut else ' ' + str(exe)}; it opens in the browser.")


if __name__ == "__main__":
    main()
