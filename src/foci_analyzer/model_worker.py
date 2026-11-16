"""The neural networks, in a process of their own.

This file runs with the Python of the models environment (`python -u model_worker.py`) and must not import
foci_analyzer: that environment only has numpy and the networks. It reads one JSON request per line from
stdin ({"op", "in", "out", "params"}; arrays travel as .npz files) and answers with one JSON line. Models stay
loaded between requests, which is what makes the preview usable.
"""

from __future__ import annotations

import json
import sys

import numpy as np

_reply = sys.stdout
sys.stdout = sys.stderr  # the libraries print progress; keep the reply channel clean
_models: dict = {}


def _device() -> str:
    import torch
    return "cuda" if torch.cuda.is_available() else "cpu"


def _version(name: str):
    try:
        from importlib.metadata import version
        return version(name)
    except Exception:  # noqa: BLE001
        return None


def ping(arrays, p):
    out = {k: _version(k) for k in ("torch", "cellpose", "stardist", "instanseg-torch", "spotiflow", "trackastra")}
    try:
        out["device"] = _device()
        if out["device"] == "cuda":
            import torch
            out["gpu"] = torch.cuda.get_device_name(0)
            out["gpu_memory"] = round(torch.cuda.get_device_properties(0).total_memory / 2**30, 1)
    except Exception:  # noqa: BLE001
        out["device"] = None
    return {}, out


def warm(arrays, p):
    """Import the libraries (the slow part of the first use) before the user asks for them."""
    import importlib

    import torch
    if torch.cuda.is_available():
        torch.zeros(1, device="cuda")
    names = {"cellpose": "cellpose.models", "instanseg": "instanseg", "stardist": "stardist.models",
             "spotiflow": "spotiflow.model", "trackastra": "trackastra.model"}
    for op in p.get("ops", []):
        try:
            importlib.import_module(names[op])
        except Exception:  # noqa: BLE001, S110  the real call reports what is missing
            pass
    return {}, {}


# ---------------------------------------------------------------- segmentation
def cellpose(arrays, p):
    from cellpose import models
    major = int((_version("cellpose") or "3").split(".")[0])
    gpu = _device() == "cuda"
    img, second = arrays["img"], arrays.get("second")
    name = p.get("model") or ("cpsam" if major >= 4 else ("nuclei" if second is None and p.get("nuclei") else "cyto3"))
    auto = major < 4 and not p.get("diameter") and name in ("nuclei", "cyto3", "cyto2", "cyto")
    key = ("cellpose", name, auto)
    if key not in _models:
        if major >= 4:
            _models[key] = models.CellposeModel(gpu=gpu) if name == "cpsam" else models.CellposeModel(gpu=gpu, pretrained_model=name)
        elif auto:
            _models[key] = models.Cellpose(gpu=gpu, model_type=name)  # with the size model: finds the diameter
        elif name in ("nuclei", "cyto3", "cyto2", "cyto"):
            _models[key] = models.CellposeModel(gpu=gpu, model_type=name)
        else:
            _models[key] = models.CellposeModel(gpu=gpu, pretrained_model=name)
    model = _models[key]
    three_d = img.shape[0] > 1
    x = img if three_d else img[0]
    kw = {"diameter": p.get("diameter") or None, "cellprob_threshold": float(p.get("cellprob", 0.0)), "batch_size": 32}
    if second is not None:
        s = second if three_d else second[0]
        x = np.stack([x, s], axis=1 if three_d else 0)  # Z, C, Y, X or C, Y, X
        if major < 4:
            kw.update(channels=[1, 2], channel_axis=1 if three_d else 0)
    elif major < 4:
        kw.update(channels=[0, 0])
    if three_d:
        kw.update(do_3D=True, z_axis=0, anisotropy=float(p.get("anisotropy", 1.0)))
    masks = model.eval(x, **kw)[0]
    return {"labels": np.asarray(masks).reshape(img.shape).astype(np.int32)}, {"model": name}


def stardist(arrays, p):
    from csbdeep.utils import normalize
    from skimage.transform import rescale, resize
    from stardist.models import StarDist2D
    name = p.get("model") or "2D_versatile_fluo"
    if ("stardist", name) not in _models:
        import os
        if os.path.isdir(name):  # a model of your own: its folder
            model = StarDist2D(None, name=os.path.basename(os.path.normpath(name)), basedir=os.path.dirname(os.path.normpath(name)))
        else:
            try:
                model = StarDist2D.from_pretrained(name)
            except OSError:  # Windows without the right to make symbolic links: use the unpacked folder itself
                base = os.path.join(os.path.expanduser("~"), ".keras", "models", "StarDist2D", name)
                model = StarDist2D(None, name=f"{name}_extracted", basedir=base)
        _models[("stardist", name)] = model
    model = _models[("stardist", name)]
    img = arrays["img"]
    scale = min(1.0, 40.0 / p["diameter"]) if p.get("diameter") else 1.0  # the model likes nuclei of about 40 px
    out = np.zeros(img.shape, np.int32)
    for z, plane in enumerate(img):
        x = normalize(rescale(plane.astype(np.float32), scale, order=1) if scale != 1 else plane, 1, 99.8)
        lab, _ = model.predict_instances(x, prob_thresh=float(p.get("prob", 0.5)))
        out[z] = resize(lab, plane.shape, order=0, preserve_range=True, anti_aliasing=False) if scale != 1 else lab
    return {"labels": out}, {"model": name}


def instanseg(arrays, p):
    from instanseg import InstanSeg
    name = p.get("model") or "fluorescence_nuclei_and_cells"
    if ("instanseg", name) not in _models:
        _models[("instanseg", name)] = InstanSeg(name, verbosity=0)
    model = _models[("instanseg", name)]
    img, second = arrays["img"], arrays.get("second")
    out = np.zeros(img.shape, np.int32)
    for z in range(img.shape[0]):
        planes = [img[z]] if second is None else [second[z], img[z]]  # nuclei first, then the cell channel
        x = np.stack(planes).astype(np.float32)
        lab, _ = model.eval_small_image(x, pixel_size=p.get("pixel_size"), target="cells" if second is not None or p.get("cells") else "nuclei",
                                        return_image_tensor=True) if max(x.shape) <= 1024 else \
            model.eval_medium_image(x, pixel_size=p.get("pixel_size"), target="cells" if second is not None or p.get("cells") else "nuclei",
                                    return_image_tensor=True)
        out[z] = np.asarray(lab.cpu() if hasattr(lab, "cpu") else lab).squeeze().astype(np.int32)
    return {"labels": out}, {"model": name}


# ---------------------------------------------------------------- foci
def spotiflow(arrays, p):
    from spotiflow.model import Spotiflow
    name = p.get("model") or "general"
    if ("spotiflow", name) not in _models:
        import os
        _models[("spotiflow", name)] = (Spotiflow.from_folder(name) if os.path.isdir(name) else Spotiflow.from_pretrained(name))
    model = _models[("spotiflow", name)]
    img = arrays["img"].astype(np.float32)
    kw = {"prob_thresh": float(p.get("floor", 0.1)), "verbose": False, "device": _device()}
    if p.get("scale") and abs(p["scale"] - 1) > 0.05:
        kw["scale"] = float(p["scale"])
    if getattr(model.config, "is_3d", False):
        pts, details = model.predict(img, **kw)
        pos = np.round(pts).astype(np.int32).reshape(-1, 3)
    else:
        flat = img.max(0)  # a 2D model: spots from the projection, their slice from the brightest value
        pts, details = model.predict(flat, **kw)
        yx = np.clip(np.round(pts).astype(np.int32).reshape(-1, 2), 0, np.array(flat.shape) - 1)
        z = img[:, yx[:, 0], yx[:, 1]].argmax(0) if len(yx) else np.zeros(0, np.int64)
        pos = np.column_stack([z, yx]).astype(np.int32)
    prob = np.asarray(details.prob, np.float32).reshape(-1)
    order = np.argsort(-prob)
    return {"pos": pos[order], "prob": prob[order]}, {"model": name}


# ---------------------------------------------------------------- tracking
def trackastra(arrays, p):
    from trackastra.model import Trackastra
    from trackastra.tracking import graph_to_ctc
    if "trackastra" not in _models:
        _models["trackastra"] = Trackastra.from_pretrained(p.get("model") or "general_2d", device=_device())
    imgs, masks = arrays["imgs"], arrays["masks"]
    graph = _models["trackastra"].track(imgs, masks, mode=p.get("mode", "greedy"))
    graph = graph[0] if isinstance(graph, tuple) else graph  # newer versions also return the masks
    table, tracked = graph_to_ctc(graph, masks, outdir=None)
    parent = {int(r.label): int(r.parent) for r in table.itertuples()}
    rows = []
    for t in range(masks.shape[0]):
        m = masks[t] > 0
        pairs = np.unique(np.stack([masks[t][m], tracked[t][m]]), axis=1).T
        rows += [(t, int(a), int(b), parent.get(int(b), 0)) for a, b in pairs if b > 0]
    return {"links": np.array(rows, np.int32).reshape(-1, 4)}, {}


OPS = {"ping": ping, "warm": warm, "cellpose": cellpose, "stardist": stardist, "instanseg": instanseg, "spotiflow": spotiflow,
       "trackastra": trackastra}


def main():
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
            arrays = dict(np.load(req["in"])) if req.get("in") else {}
            out, info = OPS[req["op"]](arrays, req.get("params") or {})
            if req.get("out"):
                np.savez(req["out"], **out)
            answer = {"ok": True, "info": info}
        except BaseException as e:
            if isinstance(e, KeyboardInterrupt):
                raise
            answer = {"ok": False, "error": f"{type(e).__name__}: {e}"}
            try:
                import torch
                torch.cuda.empty_cache()
            except Exception:  # noqa: BLE001, S110
                pass
        _reply.write(json.dumps(answer) + "\n")
        _reply.flush()


if __name__ == "__main__":
    main()
