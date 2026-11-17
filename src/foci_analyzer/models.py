"""Neural-network back-ends (Cellpose, StarDist, InstanSeg, Spotiflow, Trackastra).

The networks run in a worker process with its own Python (`model_worker.py`), so their heavy packages stay
out of the environment of the window and an existing environment can be used. Which Python:
the setting 'Python of the models', else the environment variable FOCI_ANALYZER_MODELS_PYTHON, else the
`models-env` next to the environment of Foci Analyzer, else this Python when it has PyTorch.
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import platform
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

NETWORK_PIXEL = 0.25  # µm: images sampled finer than this are made smaller for the networks

_workers: dict[str, Worker] = {}
_pythons: dict[str, str] = {}  # chosen on the Settings page: "all" and, per network, another Python
USE_GPU = True  # the Settings page can make the networks run on the processor
NETWORKS = ("cellpose", "cellpose4", "stardist", "instanseg", "spotiflow", "trackastra")
# The pretrained models offered in the window, per network: (name, description, what it needs). The empty
# name is the default of the method. Cellpose: "cellpose" runs on any version (Cellpose 4 then uses cpsam),
# "cellpose3" needs Cellpose 3, "cellpose4" Cellpose 4; they cannot be installed in one environment.
DOCS = {"stardist": "https://github.com/stardist/stardist#pretrained-models-for-2d",
        "cellpose": "https://cellpose.readthedocs.io/en/latest/models.html",
        "instanseg": "https://github.com/instanseg/instanseg",
        "spotiflow": "https://weigertlab.org/spotiflow/pretrained.html"}
CATALOG = {
    "stardist": [("", "Versatile fluorescent nuclei (2D_versatile_fluo)", "stardist")],
    "cellpose": [("", "Automatic: nuclei or cyto3 (Cellpose 3), CellposeSAM (Cellpose 4)", "cellpose"),
                 ("nuclei", "Nuclei (Cellpose 3)", "cellpose3"),
                 ("cyto3", "Cells, general (Cellpose 3)", "cellpose3"),
                 ("tissuenet_cp3", "Tissue, cells and nuclei (Cellpose 3)", "cellpose3"),
                 ("cpsam_v2", "CellposeSAM v2: any cell or nucleus (Cellpose 4)", "cellpose4"),
                 ("cpdino-vitb", "CellposeDINO ViT-B (Cellpose 4; needs the dinov3 package, see the help)", "cellpose4")],
    "instanseg": [("", "Fluorescent nuclei and cells (fluorescence_nuclei_and_cells)", "instanseg"),
                  ("brightfield_nuclei", "Nuclei in brightfield (H&E, DAB)", "instanseg")],
    "spotiflow": [("general", "General: spots from many microscopes (2D)", "spotiflow"),
                  ("hybiss", "HybISS spots (2D)", "spotiflow"),
                  ("synth_complex", "Synthetic spots with aberrations and background (2D)", "spotiflow"),
                  ("fluo_live", "Live-cell fluorescence (2D)", "spotiflow"),
                  ("synth_3d", "Synthetic 3D spots (3D stacks)", "spotiflow"),
                  ("smfish_3d", "smFISH in 3D (3D stacks)", "spotiflow")],
}


def network_of(method: str, model: str) -> str:
    """The network (environment) that runs this model: Cellpose 4 models need the Cellpose 4 environment."""
    for name, _, net in CATALOG.get(method, []):
        if name == model:
            return "cellpose" if net == "cellpose3" else net
    return method


DISTS = {"cellpose": "cellpose", "cellpose4": "cellpose", "stardist": "stardist", "instanseg": "instanseg-torch",
         "spotiflow": "spotiflow", "trackastra": "trackastra"}
_asked: dict[tuple, dict] = {}


def _versions_in(py: str) -> dict[str, str]:
    """The versions of the network packages in the environment of this Python: from its package folders (fast),
    else by asking that Python (any layout; without importing the networks, so well under a second)."""
    key = (py, Path(py).stat().st_mtime if Path(py).is_file() else 0)
    if key in _asked:
        return _asked[key]
    found: dict[str, str] = {}
    names = sorted(set(DISTS.values()) | {"torch"})
    for root in {Path(py).parent, Path(py).parent.parent}:
        for site in [root / "Lib" / "site-packages", *root.glob("lib/python3*/site-packages")]:
            for name in names:
                prefix = name.replace("-", "_")
                for d in site.glob(f"{prefix}-*.dist-info"):
                    found[name] = d.name[len(prefix) + 1:-len(".dist-info")]
    if not found and Path(py).is_file():
        code = (f"import importlib.metadata as m, json\nout = {{}}\nfor n in {names!r}:\n    try: out[n] = m.version(n)\n"
                "    except Exception: pass\nprint(json.dumps(out))")
        try:
            r = subprocess.run([py, "-c", code], check=False, env=env_for(py), capture_output=True, text=True, timeout=30,
                               creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            found = json.loads(r.stdout.strip().splitlines()[-1]) if r.returncode == 0 and r.stdout.strip() else {}
        except (OSError, ValueError, subprocess.TimeoutExpired):
            found = {}
    _asked[key] = found
    return found


def _major(version: str | None) -> int:
    try:
        return int(str(version).split(".")[0])
    except ValueError:
        return 0


def gpu_ready() -> dict[str, bool]:
    """Per network whether it will run on the graphics card: a CUDA build of PyTorch (or a Mac with Apple
    Silicon) and the graphics card switched on (Settings page). StarDist uses TensorFlow: on the processor on
    Windows and macOS."""
    mac = sys.platform == "darwin" and platform.machine() == "arm64"
    out = {}
    for op in NETWORKS:
        torch = _versions_in(python_for(op)).get("torch", "") if python_for(op) else ""
        out[op] = USE_GPU and op != "stardist" and bool(torch) and ("+cu" in torch or mac)
    return out


def installed() -> dict[str, str | None]:
    """Per network the version in the environment that would run it, or None. "cellpose" is any Cellpose (the
    automatic model), "cellpose3" Cellpose 3 (its named models), "cellpose4" Cellpose 4 (its models)."""
    out: dict[str, str | None] = {}
    for op in NETWORKS:
        py = python_for(op)
        out[op] = _versions_in(py).get(DISTS[op]) if py else None
    if _major(out["cellpose4"]) < 4:
        out["cellpose4"] = None
    out["cellpose3"] = out["cellpose"] if out["cellpose"] and _major(out["cellpose"]) < 4 else None
    return out


def catalog() -> dict:
    return {"models": {k: [{"name": n, "label": d, "network": net} for n, d, net in v] for k, v in CATALOG.items()},
            "docs": DOCS}


def set_pythons(python: str = "", per_network: dict | None = None) -> None:
    _pythons.clear()
    _pythons["all"] = str(python or "").strip().strip('"')
    for name, path in (per_network or {}).items():
        if name in NETWORKS and str(path or "").strip():
            _pythons[name] = str(path).strip().strip('"')


def python_for(op: str = "") -> str:
    """The Python that runs this network: its own, the general one, or the models environment."""
    if op in ("stardist", "cellpose4") and not _pythons.get(op):  # its own environment goes before the general one
        own = default_python(op)
        if f"{op}-env" in own:
            return own
    if op == "cellpose4" and not _pythons.get(op):  # the Cellpose of the user may be Cellpose 4 already
        py = python_for("cellpose")
        if py and _major(_versions_in(py).get("cellpose")) >= 4:
            return py
    return _pythons.get(op) or _pythons.get("all") or default_python(op)
_lock = threading.Lock()


def default_python(op: str = "") -> str:
    """The Python that holds the networks, or '' when there is none.

    StarDist has an environment of its own (`stardist-env`): it needs TensorFlow, and the TensorFlow that can
    live next to the other networks is a hundred times slower on Windows.
    """
    env = os.environ.get("FOCI_ANALYZER_MODELS_PYTHON")
    if env and Path(env).is_file():
        return env
    for folder in ([f"{op}-env"] if op in ("stardist", "cellpose4") else []) + ["models-env"]:
        for sub in ("Scripts/python.exe", "bin/python"):
            p = Path(sys.prefix).parent / folder / sub
            if p.is_file():
                return str(p)
    return sys.executable if importlib.util.find_spec("torch") else ""


def env_for(python: str) -> dict:
    """The environment variables to run this Python with. A conda environment is used without `conda activate`:
    its library folders are put on PATH as activating would (else PyTorch or numpy may fail to load a DLL)."""
    env = dict(os.environ)
    root = Path(python).parent
    if not (root / "conda-meta").is_dir():
        root = root.parent  # bin/python on macOS and Linux
    if (root / "conda-meta").is_dir():
        if sys.platform == "win32":
            dirs = [root, root / "Library" / "mingw-w64" / "bin", root / "Library" / "usr" / "bin",
                    root / "Library" / "bin", root / "Scripts", root / "bin"]
        else:
            dirs = [root / "bin"]
        env["PATH"] = os.pathsep.join([*map(str, dirs), env.get("PATH", "")])
        env["CONDA_PREFIX"] = str(root)
    return env


class Worker:
    """One running `model_worker.py`; requests go through it one at a time."""

    def __init__(self, python: str):
        self.python = python
        self.lock = threading.Lock()
        self.proc: subprocess.Popen | None = None
        self.done: set[str] = set()  # operations that ran in this process: their network is in memory
        self.dir = tempfile.mkdtemp(prefix="foci_models_")
        self.log = Path(self.dir) / "worker.log"

    def _start(self) -> None:
        script = Path(__file__).with_name("model_worker.py")
        self.proc = subprocess.Popen(
            [self.python, "-u", str(script)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=self.log.open("ab"), text=True, encoding="utf-8",
            env={**env_for(self.python), **({} if USE_GPU else {"CUDA_VISIBLE_DEVICES": "-1", "PYTORCH_ENABLE_MPS_FALLBACK": "1"})},
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))

    def call(self, op: str, arrays: dict | None = None, params: dict | None = None) -> tuple[dict, dict]:
        with self.lock:
            if self.proc is None or self.proc.poll() is not None:
                self.done.clear()
                self._start()
            req = {"op": op, "params": params or {}, "out": str(Path(self.dir) / "out.npz")}
            if arrays:
                req["in"] = str(Path(self.dir) / "in.npz")
                np.savez(req["in"], **{k: v for k, v in arrays.items() if v is not None})
            try:
                self.proc.stdin.write(json.dumps(req) + "\n")
                self.proc.stdin.flush()
                line = self.proc.stdout.readline()
            except OSError:
                line = ""
            if not line:
                tail = self.log.read_text(encoding="utf-8", errors="replace")[-600:] if self.log.is_file() else ""
                self.proc = None
                raise ValueError(f"the model process stopped ({self.python}). {tail.strip()}")
            answer = json.loads(line)
            if not answer.get("ok"):
                raise ValueError(_friendly(answer.get("error", "the model failed"), op))
            out = {}
            if Path(req["out"]).is_file():
                with np.load(req["out"]) as f:
                    out = {k: f[k] for k in f.files}
                Path(req["out"]).unlink()
            self.done.add(op)
            return out, answer.get("info", {})

    def stop(self) -> None:
        if self.proc and self.proc.poll() is None:
            self.proc.kill()


def _friendly(error: str, op: str) -> str:
    if "ModuleNotFoundError" in error or "No module named" in error:
        return f"{op} is not installed in the models environment ({error.split(':')[-1].strip()})"
    if "dinov3" in error:
        return ("CellposeDINO needs the dinov3 package of Meta in the Cellpose 4 environment, which is not installed with "
                "Foci Analyzer: pip install git+https://github.com/facebookresearch/dinov3 (in that environment)")
    if "out of memory" in error.lower():
        return f"{op}: the graphics card ran out of memory; use XY binning or 2D segmentation"
    return f"{op}: {error}"


def worker(op: str = "", python: str = "") -> Worker:
    python = python or python_for(op)
    if not python:
        raise ValueError("the neural networks are not installed: run 'install.bat' again and choose the networks, or give the Python of "
                         "an environment that has them on the Settings page")
    if not Path(python).is_file():
        raise ValueError(f"Python of the models not found: {python}")
    with _lock:
        if python not in _workers:
            _workers[python] = Worker(python)
        return _workers[python]


def warm(op: str) -> bool:
    """Whether this network ran before in a model process that is still there (so it is in memory)."""
    with _lock:
        return any(op in w.done and w.proc is not None and w.proc.poll() is None for w in _workers.values())


def stop_all() -> None:
    """End the model processes (when the program stops)."""
    with _lock:
        for w in _workers.values():
            w.stop()
        _workers.clear()


def available() -> dict:
    """What the window shows about the networks (asking the worker itself takes seconds, so not here)."""
    return {"python": default_python()}


def versions(python: str = "") -> dict:
    return worker(python=python).call("ping")[1]


STARDIST_SIZE = 40.0  # px: the nucleus diameter StarDist is given when it has to be rescaled
STARDIST_RANGE = (16, 120)  # px: nuclei in this range are left at their scale (76 px tested fine, 180 px not)


def _scaled(net: str, method: str, img: np.ndarray, second, params: dict, f: float) -> np.ndarray:
    """Labels of the network for the image scaled by `f` (below 1: smaller), back at the size of `img`."""
    if abs(f - 1) < 0.05:
        out, _ = worker(net).call(method, {"img": np.ascontiguousarray(img), "second": second}, params)
        return out["labels"]
    zoom = (1, f, f)
    small = ndi.zoom(img.astype(np.float32), zoom, order=1)
    small2 = ndi.zoom(second.astype(np.float32), zoom, order=1) if second is not None else None
    p = {**params, "diameter": params["diameter"] * f, "anisotropy": params["anisotropy"] * f,
         "pixel_size": params["pixel_size"] / f if params.get("pixel_size") else None}
    out, _ = worker(net).call(method, {"img": small, "second": small2}, p)
    lab = out["labels"]
    back = [np.minimum((np.arange(n) * f).astype(np.int64), lab.shape[i + 1] - 1) for i, n in enumerate(img.shape[1:])]
    return lab[:, back[0][:, None], back[1][None, :]]


def typical_diameter(labels: np.ndarray) -> float:
    """The median diameter (px, of a circle with the same area) of the objects, leaving out small debris."""
    area = np.bincount(labels.ravel())[1:]
    area = area[area > 0]
    if not area.size:
        return 0.0
    area = area[area >= 0.25 * np.percentile(area, 90)]
    return float(2 * np.sqrt(np.median(area) / np.pi))


def segment(img: np.ndarray, second: np.ndarray | None, cal, seg) -> np.ndarray:
    """Labels (Z, Y, X) of `img` from the chosen network; `second` is the nuclei channel next to a cell channel.

    Images sampled finer than NETWORK_PIXEL are made smaller for the network. StarDist is moreover rescaled
    so that nuclei are about STARDIST_SIZE pixels: to the Expected diameter when that is set; otherwise, when
    its nuclei come out outside STARDIST_RANGE, the best of four other scales is kept."""
    d = seg.diameter / cal.pixel_size if seg.diameter > 0 else 0
    params = {"model": seg.model, "diameter": d, "prob": seg.probability, "pixel_size": cal.pixel_size if cal.unit == "µm" else None,
              "nuclei": second is None, "anisotropy": cal.z_spacing / cal.pixel_size}
    if seg.method == "cellpose":
        p = min(max(seg.probability, 1e-3), 1 - 1e-3)
        params["cellprob"] = -math.log(p / (1 - p))  # 0.5 is the Cellpose default of 0; lower finds more
    f = cal.pixel_size / NETWORK_PIXEL if cal.unit == "µm" else 1.0
    f = f if f < 0.8 else 1.0  # only finer images are made smaller
    net = network_of(seg.method, seg.model)
    if seg.method != "stardist":
        return _scaled(net, seg.method, img, second, params, f)
    params["diameter"] = 0  # StarDist is rescaled here, not in the model process
    limit = 8000 / max(img.shape[1:])  # never larger than 8000 px
    if d > 0:
        return _scaled(net, seg.method, img, second, params, min(STARDIST_SIZE / d, limit))
    lab = _scaled(net, seg.method, img, second, params, f)
    seen = typical_diameter(lab) * f  # in pixels of the image the network saw

    def score(size):  # in the range first, then closest to STARDIST_SIZE
        return abs(size - STARDIST_SIZE) + (0 if STARDIST_RANGE[0] <= size <= STARDIST_RANGE[1] else 1e6)

    if STARDIST_RANGE[0] <= seen <= STARDIST_RANGE[1]:
        return lab
    # Too small can also mean too large: nuclei far larger than StarDist knows come out as small fragments.
    # So other scales are tried both ways, and the one with nuclei closest to STARDIST_SIZE pixels is kept.
    best = (score(seen), lab)
    for k in (0.25, 0.5, 2.0, 4.0):
        g = min(f * k, limit)
        if g * min(img.shape[1:]) < 64 or abs(g - f) < 0.05 * f:
            continue
        trial = _scaled(net, seg.method, img, second, params, g)
        best = min(best, (score(typical_diameter(trial) * g), trial), key=lambda b: b[0])
    return best[1]


def spots(raw: np.ndarray, cal, fc) -> tuple[np.ndarray, np.ndarray]:
    """Spot positions (n, 3: z, y, x) and probabilities from Spotiflow, most probable first."""
    fwhm = fc.size / cal.pixel_size
    scale = 1.0 if 1.5 <= fwhm <= 5 else min(4.0, max(0.25, 3.0 / fwhm))
    out, _ = worker("spotiflow").call("spotiflow", {"img": np.ascontiguousarray(raw)},
                           {"model": "" if fc.spotiflow_model == "general" else fc.spotiflow_model, "scale": scale, "floor": 0.1})
    return out["pos"].reshape(-1, 3), out["prob"]


def track(labels: list, images: list | None = None) -> list | None:
    """Per frame {cell: (track, parent)} from Trackastra, or None to fall back on overlap linking."""
    if images is None:
        return None
    try:
        masks = np.stack([lab.max(0) for lab in labels])
        out, _ = worker("trackastra").call("trackastra", {"imgs": np.stack(images), "masks": masks})
    except ValueError:
        return None
    frames: list[dict] = [{} for _ in labels]
    for t, cell, trk, parent in out["links"]:
        frames[int(t)][int(cell)] = (int(trk), int(parent))
    return frames
