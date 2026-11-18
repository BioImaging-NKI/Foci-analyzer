"""Cells removed by hand from the analysis.

They are kept in `removed_cells.tsv` in the results folder as points (image, frame, x, y, z in image pixels;
z = -1 for a click in a projection), not as cell numbers: after a new analysis, even with other segmentation
settings, the cell at that point is removed again. A removed cell stays in the label images but leaves the
tables, with its foci.
"""

from __future__ import annotations

import os
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from .pipeline import FrameResult

FILE = "removed_cells.tsv"
COLUMNS = ["image", "frame", "x", "y", "z"]


def load(folder: str | Path | None) -> pd.DataFrame:
    p = Path(folder or "") / FILE
    if not folder or not p.is_file():
        return pd.DataFrame(columns=COLUMNS)
    try:
        df = pd.read_csv(p, sep="\t", dtype={"image": str})
    except (OSError, ValueError):
        return pd.DataFrame(columns=COLUMNS)
    return df[COLUMNS] if set(COLUMNS) <= set(df.columns) else pd.DataFrame(columns=COLUMNS)


def save(folder: str | Path, df: pd.DataFrame) -> None:
    p = Path(folder) / FILE
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + ".tmp")
    df[COLUMNS].to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, p)


def points(df: pd.DataFrame, image_id: str, t: int) -> list[tuple[float, float, int]]:
    rows = df[(df["image"] == image_id) & (df["frame"] == t + 1)]
    return [(float(r.x), float(r.y), int(r.z)) for r in rows.itertuples()]


def label_at(r: FrameResult, x: float, y: float, z: int) -> int:
    """The cell under a point given in image pixels (0 = none)."""
    lab = r.labels
    xi = int((x - r.offset[0] - r.border) // r.binning)
    yi = int((y - r.offset[1] - r.border) // r.binning)
    if not (0 <= yi < lab.shape[1] and 0 <= xi < lab.shape[2]):
        return 0
    if lab.shape[0] == 1:
        return int(lab[0, yi, xi])
    if 0 <= z < lab.shape[0]:
        return int(lab[z, yi, xi])
    column = lab[:, yi, xi]  # clicked in a projection: the cell that is there in any slice
    hit = column[column > 0]
    return int(np.bincount(hit).argmax()) if hit.size else 0


def apply(r: FrameResult, pts: list[tuple[float, float, int]]) -> FrameResult:
    """A copy of the result without the removed cells (and their foci) in the tables."""
    ids = sorted({label_at(r, *p) for p in pts} - {0})
    if not ids:
        return replace(r, removed=[])
    cells = r.cells[~r.cells["cell"].isin(ids)] if "cell" in r.cells else r.cells
    foci = r.foci[~r.foci["cell"].isin(ids)] if "cell" in r.foci else r.foci
    return replace(r, cells=cells.reset_index(drop=True), foci=foci.reset_index(drop=True), removed=ids)


def set_removed(df: pd.DataFrame, r: FrameResult, image_id: str, t: int, pts: list, remove: bool) -> tuple[pd.DataFrame, int]:
    """Remove the cells under the points (x, y, z) from the analysis, or put them back. Returns the table and how
    many cells changed."""
    changed = 0
    for x, y, z in pts:
        cell = label_at(r, float(x), float(y), int(z))
        if cell == 0:
            continue
        here = (df["image"] == image_id) & (df["frame"] == t + 1)
        same = [i for i, row in df[here].iterrows() if label_at(r, row.x, row.y, int(row.z)) == cell]
        if remove and not same:
            row = pd.DataFrame([[image_id, t + 1, round(float(x), 1), round(float(y), 1), int(z)]], columns=COLUMNS)
            df = row if df.empty else pd.concat([df, row], ignore_index=True)
            changed += 1
        elif not remove and same:
            df = df.drop(index=same).reset_index(drop=True)
            changed += 1
    return df, changed


def toggle(df: pd.DataFrame, r: FrameResult, image_id: str, t: int, x: float, y: float, z: int) -> tuple[pd.DataFrame, bool]:
    """Remove the cell under (x, y, z), or put it back when it was removed. Returns the table and whether it is
    now removed."""
    cell = label_at(r, x, y, z)
    if cell == 0:
        raise ValueError("there is no cell at this place")
    here = (df["image"] == image_id) & (df["frame"] == t + 1)
    same = [i for i, row in df[here].iterrows() if label_at(r, row.x, row.y, int(row.z)) == cell]
    if same:
        return df.drop(index=same).reset_index(drop=True), False
    row = pd.DataFrame([[image_id, t + 1, round(x, 1), round(y, 1), z]], columns=COLUMNS)
    return (row if df.empty else pd.concat([df, row], ignore_index=True)), True
