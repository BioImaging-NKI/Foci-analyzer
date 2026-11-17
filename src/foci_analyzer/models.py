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
# The pretrained models offered in the window, per network: (name, description, network that runs it).
# Cellpose 3 and 4 need environments of their own; the empty name is the default of the method.
DOCS = {"stardist": "https://github.com/stardist/stardist#pretrained-models-for-2d",
        "cellpose": "https://cellpose.readthedocs.io/en/latest/models.html",
        "instanseg": "https://github.com/instanseg/instanseg",
        "spotiflow": "https://weigertlab.org/spotiflow/pretrained.html"}
CATALOG = {
    "stardist": [("", "Versatile fluorescent nuclei (2D_versatile_fluo)", "stardist"),
                 ("2D_paper_dsb2018", "Nuclei of the Data Science Bowl 2018", "stardist")],
    "cellpose": [("", "Automatic: nuclei, or cyto3 with a cell channel (Cellpose 3)", "cellpose"),
                 ("nuclei", "Nuclei (Cellpose 3)", "cellpose"),
                 ("cyto3", "Cells, general (Cellpose 3)", "cellpose"),
                 ("tissuenet_cp3", "Tissue, cells and nuclei (Cellpose 3)", "cellpose"),
                 ("cpsam_v2", "CellposeSAM v2: any cell or nucleus (Cellpose 4, slow without GPU)", "cellpose4"),
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
            return net
    return method


def installed() -> dict[str, str | None]:
    """Per network the version in the environment that would run it, or None; from the package folders, so
    without starting that Python (which takes seconds)."""
    dist = {"cellpose": "cellpose", "cellpose4": "cellpose", "stardist": "stardist", "instanseg": "instanseg_torch",
            "spotiflow": "spotiflow", "trackastra": "trackastra"}
    out: dict[str, str | None] = {}
    for op in NETWORKS:
        py, version = python_for(op), None
        if py and Path(py).is_file():
            for root in {Path(py).parent, Path(py).parent.parent}:
                for site in [root / "Lib" / "site-packages", *root.glob("lib/python3*/site-packages")]:
                    for d in site.glob(f"{dist[op]}-*.dist-info"):
                        version = d.name[len(dist[op]) + 1:-len(".dist-info")]
        if version and op.startswith("cellpose") and (version.split(".")[0] >= "4") != (op == "cellpose4"):
            version = None  # Cellpose 3 models need Cellpose 3, the Cellpose 4 models Cellpose 4
        out[op] = version
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
            env={**os.environ, **({} if USE_GPU else {"CUDA_VISIBLE_DEVICES": "-1", "PYTORCH_ENABLE_MPS_FALLBACK": "1"})},
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


def segment(img: np.ndarray, second: np.ndarray | None, cal, seg) -> np.ndarray:
    """Labels (Z, Y, X) of `img` from the chosen network; `second` is the nuclei channel next to a cell channel."""
    d = seg.diameter / cal.pixel_size if seg.diameter > 0 else 0
    params = {"model": seg.model, "diameter": d, "prob": seg.probability, "pixel_size": cal.pixel_size if cal.unit == "µm" else None,
              "nuclei": second is None, "anisotropy": cal.z_spacing / cal.pixel_size}
    if seg.method == "cellpose":
        p = min(max(seg.probability, 1e-3), 1 - 1e-3)
        params["cellprob"] = -math.log(p / (1 - p))  # 0.5 is the Cellpose default of 0; lower finds more
    f = cal.pixel_size / NETWORK_PIXEL if cal.unit == "µm" else 1.0
    net = network_of(seg.method, seg.model)
    if f < 0.8:  # finer than the networks need: give them a smaller image and bring the labels back
        zoom = (1, f, f)
        small = ndi.zoom(img.astype(np.float32), zoom, order=1)
        small2 = ndi.zoom(second.astype(np.float32), zoom, order=1) if second is not None else None
        params.update(diameter=d * f, pixel_size=cal.pixel_size / f, anisotropy=params["anisotropy"] * f)
        out, _ = worker(net).call(seg.method, {"img": small, "second": small2}, params)
        lab = out["labels"]
        back = [np.minimum((np.arange(n) * f).astype(np.int64), lab.shape[i + 1] - 1) for i, n in enumerate(img.shape[1:])]
        return lab[:, back[0][:, None], back[1][None, :]]
    out, _ = worker(net).call(seg.method, {"img": np.ascontiguousarray(img), "second": second}, params)
    return out["labels"]


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
