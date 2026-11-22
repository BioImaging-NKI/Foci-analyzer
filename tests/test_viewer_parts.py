"""Segmenting only, reusing a segmentation for a part of the image, and the outlines of the cells."""

import threading
from collections import OrderedDict

import numpy as np

from foci_analyzer import analysis as an
from foci_analyzer import export, pipeline
from foci_analyzer.server import App
from foci_analyzer.stages import StageCache
from test_analysis import settings, synthetic


def test_segment_only_then_estimate():
    data, info, _ = synthetic()
    s = settings()
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s, segment_only=True)
    assert r.channels == [] and r.foci.empty and len(r.cells) == 2 and r.labels.max() == 2
    prepared, cal, _ = pipeline._prepare(data, info, s)
    est = an.estimate_foci(prepared[1], r.labels, cal, s.foci_general, s.foci[0])
    assert est["size"] > 0 and est["foci"] >= 6


def test_given_labels_are_used():
    data, info, _ = synthetic()
    s = settings()
    known = np.zeros((1, 128, 160), np.int32)
    known[0, 20:60, 20:70] = 7
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s, labels=known)
    assert r.labels.max() == 1 and (r.labels > 0).sum() == 40 * 50  # numbered again, the rest as found


def fake_app():
    app = object.__new__(App)
    app.lock, app._segs = threading.Lock(), OrderedDict()
    return app


def test_part_inside_a_segmented_part_is_cut_out():
    data, info, _ = synthetic()
    s = settings()
    s.image.xy_binning = 2
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s)
    app = fake_app()
    key = ("img", 0, "settings")
    app._segs[key + (0, 0)] = (r.labels, r.border, r.binning)
    lab, (x0, y0) = app._known_cells(key, 11, 21, 151, 121, info, s)  # an odd corner: moved to an even one
    assert (x0, y0) == (10, 20)
    assert lab.shape[1:] == ((121 - 20) // 2, (151 - 10) // 2)
    assert np.array_equal(lab, r.labels[:, 10:10 + lab.shape[1], 5:5 + lab.shape[2]])
    none, corner = app._known_cells(key, 100, 100, 400, 300, info, s)  # not inside: segment again
    assert none is None and corner == (100, 100)
    assert app._known_cells(("img", 0, "other settings"), 11, 21, 151, 121, info, s)[0] is None


def test_every_piece_of_a_cell_has_an_outline():
    lab = np.zeros((40, 40), np.int32)
    lab[5:15, 5:15] = 3
    lab[25:35, 20:30] = 3  # the same cell, in two pieces
    lab[18:22, 2:6] = 4
    found = [i for i, _ in export._contours(lab)]
    assert sorted(found) == [3, 3, 4]


def test_spline_rounds_off_the_pixel_steps_but_keeps_the_shape():
    yy, xx = np.mgrid[:120, :120]
    lab = ((yy - 59.5) ** 2 + (xx - 59.5) ** 2 < 48**2).astype(np.int32)
    lab[50:70, 100:] = 0  # a dent, which the outline must follow
    (_, pts), = list(export._contours(lab))
    smooth = export.spline_outline(pts)
    from scipy.spatial import cKDTree
    assert cKDTree(pts).query(smooth)[0].max() < 2.0  # no corner cut off

    def turning(p):  # total turning angle in full turns: 1 for a smooth convex outline, more for steps and dents
        d = np.diff(np.vstack([p, p[:1]]), axis=0)
        a = np.arctan2(d[:, 1], d[:, 0])
        return np.abs((np.diff(np.r_[a, a[:1]]) + np.pi) % (2 * np.pi) - np.pi).sum() / (2 * np.pi)

    assert turning(pts) > 5 and turning(smooth) < 2.3  # the dent adds one turn of its own corners

    def area(p):
        return 0.5 * abs(np.dot(p[:, 0], np.roll(p[:, 1], 1)) - np.dot(p[:, 1], np.roll(p[:, 0], 1)))

    assert abs(area(smooth) / area(pts) - 1) < 0.03
