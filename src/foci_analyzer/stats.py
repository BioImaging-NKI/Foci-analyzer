"""Comparing conditions (plan.md, section 8).

The unit of a test is a cell, an image or a folder (one folder = one experiment). Cells of one image are not
independent, so the tests on image or folder means are the careful ones. Rank tests are used per cell
(Kruskal-Wallis, Mann-Whitney), Welch tests on means; pairwise p-values get the Holm correction.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats as st


def holm(p: list[float]) -> list[float]:
    """Holm-corrected p-values, in the order given."""
    order = np.argsort(p)
    out, worst = [0.0] * len(p), 0.0
    for rank, i in enumerate(order):
        worst = max(worst, min(1.0, (len(p) - rank) * p[i]))
        out[i] = worst
    return out


def groups(df: pd.DataFrame, measure: str, unit: str) -> dict[str, np.ndarray]:
    """The values per condition, one per cell or one mean per image or folder."""
    df = df[["condition", "image", "folder", measure]].dropna(subset=[measure])
    if unit == "image":
        df = df.groupby(["condition", "folder", "image"], sort=False, as_index=False)[measure].mean()
    elif unit == "folder":
        df = df.groupby(["condition", "folder"], sort=False, as_index=False)[measure].mean()
    return {c: g[measure].to_numpy(float) for c, g in df.groupby("condition", sort=False)}


def compare(df: pd.DataFrame, measure: str, unit: str = "cell") -> dict:
    if measure not in df.columns or not pd.api.types.is_numeric_dtype(df[measure]):
        raise ValueError(f"'{measure}' is not a number column")
    g = groups(df, measure, unit)
    summary = [{"condition": c, "n": int(v.size), "mean": float(v.mean()) if v.size else None,
                "sd": float(v.std(ddof=1)) if v.size > 1 else None,
                "median": float(np.median(v)) if v.size else None} for c, v in g.items()]
    usable = {c: v for c, v in g.items() if v.size >= 2 and np.ptp(v) >= 0}
    out = {"measure": measure, "unit": unit, "summary": summary, "overall": None, "pairs": [], "note": ""}
    if len(usable) < 2:
        out["note"] = "A test needs two conditions with at least two values each."
        return out
    rank = unit == "cell"
    out["test"] = "Mann-Whitney" if rank else "Welch t-test"
    if len(usable) > 2:
        try:
            res = st.kruskal(*usable.values()) if rank else st.f_oneway(*usable.values())
            out["overall"] = {"test": "Kruskal-Wallis" if rank else "one-way ANOVA", "p": float(res.pvalue)}
        except ValueError:
            pass
    pairs = []
    for a, b in combinations(usable, 2):
        try:
            if rank:
                p = st.mannwhitneyu(usable[a], usable[b], alternative="two-sided").pvalue
            else:
                p = st.ttest_ind(usable[a], usable[b], equal_var=False).pvalue
        except ValueError:
            p = np.nan
        pairs.append([a, b, float(p)])
    good = [i for i, x in enumerate(pairs) if np.isfinite(x[2])]
    adjusted = holm([pairs[i][2] for i in good])
    for i, x in enumerate(pairs):
        out["pairs"].append({"a": x[0], "b": x[1], "p": x[2] if np.isfinite(x[2]) else None,
                             "p_holm": adjusted[good.index(i)] if i in good else None})
    if not rank and min(v.size for v in usable.values()) < 3:
        out["note"] = "Fewer than three values in a condition: the test has little power."
    return out
