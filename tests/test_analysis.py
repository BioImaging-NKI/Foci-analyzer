"""Tests of the analysis on a made-up image with known nuclei and foci."""

import numpy as np
import pandas as pd
import pytest
import tifffile

from foci_analyzer import analysis as an
from foci_analyzer import batch, config, pipeline, stats, tracking
from foci_analyzer.experiment import ImageEntry
from foci_analyzer.io import ImageInfo
from foci_analyzer.stages import StageCache

PX, ZS = 0.2, 0.6


def synthetic(seed=0, slices=5, shift=0.0):
    """Two nuclei (channel 1) with 6 and 3 foci in channel 2; channel 3 has the same foci, moved by `shift` px."""
    rng = np.random.default_rng(seed)
    z, y, x = np.mgrid[:slices, :128, :160].astype(np.float32)
    nuclei = np.zeros((slices, 128, 160), np.float32)
    for cy, cx in ((40, 45), (85, 110)):
        nuclei += 800 * (((y - cy) / 26) ** 2 + ((x - cx) / 30) ** 2 < 1)
    spots = [(2, 30, 35), (2, 45, 60), (1, 52, 38), (3, 28, 52), (2, 40, 30), (2, 55, 55), (2, 80, 100), (1, 95, 120), (3, 90, 105)]

    def foci(dx):
        out = np.zeros_like(nuclei)
        for sz, sy, sx in spots:
            out += 600 * np.exp(-(((y - sy) ** 2 + (x - sx - dx) ** 2) / (2 * 1.2**2) + (z - (sz if slices > 1 else 0)) ** 2 / (2 * 0.8**2)))
        return out

    data = np.stack([nuclei, 60 * (nuclei > 0) + foci(0), 60 * (nuclei > 0) + foci(shift)])
    data = rng.poisson(data + 20).astype(np.uint16)
    info = ImageInfo("x.tif", "x", 1, 3, slices, 128, 160, "uint16", PX, ZS if slices > 1 else None, "µm")
    return data, info, spots


def settings(**foci):
    s = config.Settings()
    s.segmentation.exclude_on_edges, s.segmentation.method = False, "threshold"
    s.foci = [config.FociChannelSettings(channel=2, **foci), config.FociChannelSettings(channel=3, **foci)]
    return s


@pytest.mark.parametrize("dimension", ["3D", "max projection", "extended depth of focus", "sum projection"])
@pytest.mark.parametrize("flt", ["DoG", "LoG", "wavelet"])
def test_finds_the_foci(dimension, flt):
    data, info, spots = synthetic()
    s = settings(filter=flt)
    s.foci_general.dimension = dimension
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s)
    assert len(r.cells) == 2
    assert sorted(r.cells["foci_ch2"]) == [3, 6]
    f = r.foci[r.foci.channel == 2]
    for _, sy, sx in spots:  # every focus within a pixel of where it was put
        assert np.min(np.hypot(f.x / PX - sx, f.y / PX - sy)) < 1.0
    assert (r.cells["fraction_ch2_with_ch3"] == 1).all()
    assert (r.cells["pearson_ch2_ch3"] > 0.5).all()


def test_colocalization_ends_when_foci_are_apart():
    data, info, _ = synthetic(shift=9.0)
    s = settings()
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s)
    assert (r.cells["fraction_ch2_with_ch3"] == 0).all()
    s.colocalization.method, s.colocalization.max_distance = "centroid distance", 2.5
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s)
    assert (r.cells["fraction_ch2_with_ch3"] == 1).all()  # 9 px = 1.8 µm is inside 2.5 µm
    s.colocalization.chance_repeats = 20
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s)
    assert r.cells["chance_ch2_with_ch3"].between(0, 1).all()


def test_threshold_change_repeats_only_the_last_stages():
    data, info, _ = synthetic()
    s, cache = settings(), StageCache()
    pipeline.analyze_frame(cache, "k", data, info, s)
    cache.computed.clear()
    s.foci[0].sensitivity = 12.0
    pipeline.analyze_frame(cache, "k", data, info, s)
    assert cache.computed == ["extents", "tables"]
    cache.computed.clear()
    s.foci[1].size = 0.7
    pipeline.analyze_frame(cache, "k", data, info, s)
    assert "segment" not in cache.computed and cache.computed.count("filter") == 1


def test_half_maximum_extent_does_not_depend_on_the_threshold():
    data, info, _ = synthetic()
    sizes = []
    for thr in (6.0, 15.0):
        r = pipeline.analyze_frame(StageCache(), "k", data, info, settings(sensitivity=thr))
        sizes.append(r.foci[r.foci.channel == 2]["size"].sum())
    assert sizes[0] == pytest.approx(sizes[1], rel=0.05)
    t = [pipeline.analyze_frame(StageCache(), "k", data, info, settings(sensitivity=thr, extent="threshold"))
         .foci["size"].sum() for thr in (6.0, 15.0)]
    assert t[0] > 1.2 * t[1]


def test_scopes_binning_and_border():
    data, info, _ = synthetic()
    for scope in ("per cell", "fixed value"):
        r = pipeline.analyze_frame(StageCache(), "k", data, info, settings(threshold_scope=scope))
        assert r.cells["foci_ch2"].sum() in (9, 10)  # a fixed default of 5 x noise may let one more through
    s = settings(size=0.6)
    s.image.xy_binning, s.image.crop_border = 2, 4
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s)
    assert r.cal.pixel_size == pytest.approx(2 * PX) and r.labels.shape[-1] == (160 - 8) // 2
    assert r.cells["foci_ch2"].sum() >= 7


def test_no_cells_and_labels_from_file(tmp_path):
    data, info, _ = synthetic(slices=1)
    s = settings()
    s.image.nuclei_channel = 0
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s)
    assert len(r.cells) == 1 and r.cells["foci_ch2"][0] >= 9 and r.cells["check"][0] == ""
    lab = np.zeros((128, 160), np.uint16)
    lab[10:70, 10:80] = 7
    tifffile.imwrite(tmp_path / "x.tif", lab)
    s = settings()
    s.segmentation.method, s.segmentation.labels_folder = "labels from file", str(tmp_path)
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s)
    assert len(r.cells) == 1 and r.cells["foci_ch2"][0] >= 6


def test_expand_and_filter_labels():
    cal = an.Calibration(PX, ZS, "µm")
    lab = np.zeros((1, 60, 60), np.int32)
    lab[0, 20:30, 20:30] = 1
    lab[0, 0:4, 0:4] = 2
    grown = an.expand(lab, cal, 1.0)
    assert grown[0, 17, 25] == 1 and grown[0, 10, 25] == 0
    assert (an.expand(lab, cal, -1) > 0).all()
    kept = an.filter_labels(lab, cal, 1.5, 50, True)
    assert kept.max() == 1 and kept[0, 25, 25] == 1 and kept[0, 1, 1] == 0


def test_holm_and_compare():
    assert stats.holm([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.06, 0.06])
    rng = np.random.default_rng(0)
    rows = [{"condition": c, "folder": f"{c}{e}", "image": f"im{i}", "foci": v}
            for c, mu in (("a", 5), ("b", 9), ("c", 5)) for e in range(3) for i in range(2)
            for v in rng.poisson(mu, 15)]
    df = pd.DataFrame(rows)
    out = stats.compare(df, "foci", "folder")
    assert [s["n"] for s in out["summary"]] == [3, 3, 3] and out["test"] == "Welch t-test"
    p = {(x["a"], x["b"]): x["p_holm"] for x in out["pairs"]}
    assert p[("a", "b")] < 0.05 and p[("a", "c")] > 0.05
    assert stats.compare(df, "foci", "cell")["summary"][0]["n"] == 90
    with pytest.raises(ValueError, match="not a number"):
        stats.compare(df, "image", "cell")


def test_overlap_tracking():
    a = np.zeros((1, 40, 40), np.int32)
    a[0, 5:15, 5:15], a[0, 25:35, 25:35] = 1, 2
    b = np.zeros_like(a)
    b[0, 26:36, 26:36], b[0, 6:11, 6:16], b[0, 11:16, 6:16] = 1, 2, 3  # moved; the other one split in two
    links = tracking.link_overlap([a, b])
    assert links[0] == {1: (1, 0), 2: (2, 0)}
    assert links[1][1] == (2, 0)
    assert {links[1][2][0], links[1][3][0]} == {1, 3} and {links[1][2][1], links[1][3][1]} == {0, 1}


def test_batch_writes_tables(tmp_path):
    data, _, _ = synthetic()
    for name in ("one", "two"):
        tifffile.imwrite(tmp_path / f"{name}.tif", data.transpose(1, 0, 2, 3), imagej=True, metadata={
            "axes": "ZCYX", "spacing": ZS, "unit": "micron"}, resolution=(1 / PX, 1 / PX))
    images = [ImageEntry(str(tmp_path / f"{n}.tif"), n, "cond", str(tmp_path), 0, str(tmp_path / f"{n}.tif"))
              for n in ("one", "two")]
    import threading
    results, progress = batch.Results(), batch.Progress()
    batch.run(images, settings(), str(tmp_path / "out"), results, progress, threading.Lock())
    assert not progress.errors and progress.done == 2
    assert len(results.cells) == 4 and results.cells["foci_ch2"].sum() == 18
    again = batch.Results()
    assert again.load(tmp_path / "out") and len(again.foci) == len(results.foci)
    assert again.meta["settings"] == config.to_dict(settings())


def test_preview_overlays_parts_and_cancel(tmp_path, monkeypatch):
    """The window's side of an analysis: preview, overlays of every kind, a part of the image, cancel."""
    monkeypatch.setenv("FOCI_ANALYZER_HOME", str(tmp_path / "home"))
    from foci_analyzer.server import ApiError, App
    data, _, _ = synthetic()
    path = tmp_path / "one.tif"
    tifffile.imwrite(path, data.transpose(1, 0, 2, 3), imagej=True, metadata={"axes": "ZCYX", "spacing": ZS, "unit": "micron"},
                     resolution=(1 / PX, 1 / PX))
    app = App()
    app.set_conditions([{"name": "c", "folders": [str(path)], "pattern": "*", "recursive": False}])
    image = app.images[0].id
    with pytest.raises(ApiError):  # an overlay never starts an analysis itself
        app.overlay({"path": image, "kind": "cells"})
    s = config.to_dict(settings())
    s["colocalization"]["enabled"] = True
    p = app.preview({"path": image, "settings": s})
    assert p["cells"] == 2 and [c["foci"] for c in p["channels"]] == [9, 9] and p["part"] is None
    assert {(c["a"], c["b"]): c["with"] for c in p["coloc"]} == {(2, 3): 9, (3, 2): 9}
    assert all(c["chance"] is not None and c["chance"] < 0.5 for c in p["coloc"])
    for kind, extra in (("cells", {}), ("foci", {"ch": "2"}), ("coloc", {}), ("overlap", {}), ("value", {"col": "foci_ch2"})):
        raw, w, h = app.overlay({"path": image, "kind": kind, "z": "mip", **extra})
        plane = np.frombuffer(raw, "<u2").reshape(h, w)
        assert (w, h) == (160, 128) and plane.any(), kind
        if kind == "value":  # the cell with 6 foci is the top of the colour scale, the one with 3 the bottom
            assert set(np.unique(plane)) == {0, 1, 255}
    cand = np.frombuffer(app.candidates({"path": image, "ch": "2"}), "<f4").reshape(-1, 6)
    assert (cand[:, 5] == 2).sum() == 9  # every focus colocalizes
    # only the left half: one nucleus, and everything is placed where it is in the whole image
    part = app.preview({"path": image, "settings": s, "box": [0, 0, 80, 128]})
    assert part["cells"] == 1 and part["offset"] == [0, 0] and part["part"][2:] == [80, 128]
    right = app.preview({"path": image, "settings": s, "box": [72, 0, 160, 128]})
    assert right["cells"] == 1 and right["offset"] == [72, 0]
    cand = np.frombuffer(app.candidates({"path": image, "ch": "2"}), "<f4").reshape(-1, 6)
    assert cand[cand[:, 5] > 0, 0].min() > 90  # x of the foci of the right nucleus, in image pixels
    assert app.cancel() == {"stopped": True}
    s["foci"][0]["size"] = 0.7  # something new to compute: it is stopped before it starts
    app.set_settings(s)
    app.cache.cancelled = True
    with pytest.raises(ApiError, match="stopped"):
        app.frame(image, 0)
    assert app.preview({"path": image, "settings": s})["cells"] == 2  # a new request clears the stop


def test_wavelet_counts_like_dog():
    data, info, _ = synthetic(seed=3)
    counts = {}
    for flt in ("DoG", "wavelet"):
        r = pipeline.analyze_frame(StageCache(), "w", data, info, settings(filter=flt))
        counts[flt] = [len(c.kept) for c in r.channels]
    assert counts["wavelet"] == counts["DoG"] == [9, 9]
