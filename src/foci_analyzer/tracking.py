"""Following cells over the frames of a time-lapse (plan.md, section 5.7).

`overlap`: a cell continues the track of the cell in the previous frame that it overlaps most. When two
cells claim the same earlier cell (a division, or a split in the segmentation) the larger overlap keeps the
track and the other starts a new track with `parent` set.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def link_overlap(labels: list[np.ndarray]) -> list[dict[int, tuple[int, int]]]:
    """Per frame: cell number → (track, parent track or 0)."""
    out, next_track, prev, prev_map = [], 1, None, {}
    for lab in labels:
        foot = lab.max(0)
        n = int(foot.max())
        current: dict[int, tuple[int, int]] = {}
        best: dict[int, tuple[int, int]] = {}  # earlier cell → (overlap, cell) of its best follower
        claims: dict[int, tuple[int, int]] = {}  # cell → (earlier cell, overlap)
        if prev is not None and n and prev.shape == foot.shape:
            m = (foot > 0) & (prev > 0)
            pairs, count = np.unique(foot[m].astype(np.int64) * (int(prev.max()) + 1) + prev[m], return_counts=True)
            for pair, k in zip(pairs, count, strict=True):
                cell, before = divmod(int(pair), int(prev.max()) + 1)
                if k > claims.get(cell, (0, 0))[1]:
                    claims[cell] = (before, int(k))
            for cell, (before, k) in claims.items():
                if k > best.get(before, (0, 0))[0]:
                    best[before] = (k, cell)
        for cell in range(1, n + 1):
            before = claims.get(cell, (0, 0))[0]
            if before and best[before][1] == cell:
                current[cell] = (prev_map[before][0], prev_map[before][1])
            else:
                current[cell] = (next_track, prev_map[before][0] if before else 0)
                next_track += 1
        out.append(current)
        prev, prev_map = foot, current
    return out


def add_tracks(cells: pd.DataFrame, foci: pd.DataFrame, labels: list[np.ndarray], settings,
               images: list[np.ndarray] | None = None) -> str:
    """Add the columns `track` and `parent` to the cell table and `track` to the foci table.

    Returns a remark when Trackastra was asked for but could not be used.
    """
    links, remark = None, ""
    if settings.tracking == "trackastra":
        from . import models
        links = models.track(labels, images)
        if links is None:
            remark = "Trackastra could not be used; cells were followed by overlap"
    if links is None:
        links = link_overlap(labels)
    track = {(t + 1, cell): v for t, frame in enumerate(links) for cell, v in frame.items()}
    keys = list(zip(cells["frame"], cells["cell"], strict=True))
    cells["track"] = [track.get(k, (0, 0))[0] for k in keys]
    cells["parent"] = [track.get(k, (0, 0))[1] for k in keys]
    if len(foci):
        foci["track"] = [track.get(k, (0, 0))[0] for k in zip(foci["frame"], foci["cell"], strict=True)]
    return remark
