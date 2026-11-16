"""Check an installation: every installed network segments a made-up image of nuclei, and Spotiflow finds spots.

    <program folder>/.venv/Scripts/python.exe installer/check.py      (Windows)
    <program folder>/.venv/bin/python installer/check.py              (macOS, Linux)

Prints per network where it runs (graphics card or processor) and whether it worked. Used by the automatic
builds on GitHub, and useful for anyone who wants to know whether an installation works.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from foci_analyzer import analysis as an
from foci_analyzer import config, models


def nuclei_image() -> tuple[np.ndarray, list]:
    """A 2D image (Z = 1) with 12 bright round nuclei of about 10 µm and some foci, 0.25 µm per pixel."""
    rng = np.random.default_rng(1)
    y, x = np.mgrid[:512, :512]
    img = np.zeros((512, 512), np.float32)
    centres = [(64 + 128 * (i % 4), 85 + 170 * (i // 4)) for i in range(12)]
    for cx, cy in centres:
        img += 900 * (((x - cx) / 22) ** 2 + ((y - cy) / 19) ** 2 < 1)
    spots = img * 0
    for cx, cy in centres:
        for _ in range(5):
            sx, sy = cx + rng.integers(-12, 12), cy + rng.integers(-10, 10)
            spots += 1500 * np.exp(-((x - sx) ** 2 + (y - sy) ** 2) / (2 * 1.3**2))
    img = rng.poisson(img + 40).astype(np.float32)
    return img[None], [rng.poisson(spots + img * 0.2 + 30).astype(np.float32)[None]]


def main() -> int:
    img, (spots,) = nuclei_image()
    cal = an.Calibration(0.25, 0.25, "µm")
    failed = 0
    for name in models.NETWORKS:
        python = models.python_for(name)
        if not python:
            print(f"{name:10s}  not installed")
            continue
        try:
            info = models.versions(python)
            version = info.get("instanseg-torch" if name == "instanseg" else name)
            if not version:
                print(f"{name:10s}  not installed in {python}")
                continue
            where = {"cuda": f"graphics card ({info.get('gpu')})", "mps": "Apple graphics chip",
                     "cpu": "processor"}.get(info.get("device"), "?")
            if name == "stardist":
                where = "processor (TensorFlow)"
            t = time.perf_counter()
            if name == "spotiflow":
                fc = config.FociChannelSettings(method="spotiflow", size=0.5)
                pos, _ = models.spots(spots, cal, fc)
                found, ok = f"{len(pos)} spots", 30 <= len(pos) <= 120
            elif name == "trackastra":
                lab = an.segment_threshold(img, cal, 4.0)
                frames = models.track([lab, lab], [img[0], img[0]])
                found, ok = "tracked" if frames else "no tracks", bool(frames) and len(frames[1]) >= 10
            else:
                seg = config.SegmentationSettings(method=name)
                lab = models.segment(img, None, cal, seg)
                n = len(np.unique(lab)) - 1
                found, ok = f"{n} nuclei", 10 <= n <= 14
            print(f"{name:10s}  {version:10s} {where:40s} {found:12s} {time.perf_counter() - t:5.1f} s  {'ok' if ok else 'WRONG'}")
            failed += not ok
        except ValueError as e:
            print(f"{name:10s}  FAILED: {e}")
            failed += 1
    models.stop_all()
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
