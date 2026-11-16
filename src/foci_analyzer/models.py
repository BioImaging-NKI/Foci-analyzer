"""Neural-network back-ends (Cellpose, StarDist, InstanSeg, Spotiflow). Filled in with milestones 3 and 4."""

from __future__ import annotations

import numpy as np


def segment(img: np.ndarray, second: np.ndarray | None, cal, seg) -> np.ndarray:
    raise ValueError(f"the segmentation method '{seg.method}' is not available yet; use 'threshold' or "
                     "'labels from file'")


def spots(raw: np.ndarray, cal, fc) -> tuple[np.ndarray, np.ndarray]:
    raise ValueError("Spotiflow is not available yet; use 'watershed'")


def available() -> dict:
    """Which back-ends can be used in this installation."""
    return {}


def track(labels: list) -> list | None:
    """Tracks from Trackastra, or None to fall back on overlap linking."""
    return None
