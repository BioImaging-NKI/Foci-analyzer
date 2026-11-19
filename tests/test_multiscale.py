"""The multi-scale foci method: foci of very different sizes, each found once and at about its own size."""

import itertools

import numpy as np
import pytest

from foci_analyzer import analysis as an
from foci_analyzer import config, pipeline
from foci_analyzer.io import ImageInfo
from foci_analyzer.stages import StageCache
from test_analysis import PX as PX0
from test_analysis import synthetic

PX = 0.1
SIZES = (0.3, 0.6, 1.2, 2.0)  # µm, full width at half maximum


def mixed(seed=0):
    """One large nucleus (channel 1) with three foci of every size in SIZES (channel 2), on a grid."""
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[:300, :300].astype(np.float32)
    nucleus = 800.0 * (((y - 150) / 140) ** 2 + ((x - 150) / 140) ** 2 < 1)
    foci, truth = np.zeros_like(nucleus), []
    spots = [(70 + 55 * i, 70 + 55 * j) for i in range(3) for j in range(4)]
    for (sy, sx), size in zip(spots, [s for _ in range(3) for s in SIZES], strict=True):
        sigma = size / PX / an.FWHM
        foci += 500 * np.exp(-((y - sy) ** 2 + (x - sx) ** 2) / (2 * sigma**2))
        truth.append((sy, sx, size))
    data = np.stack([nucleus, 80 * (nucleus > 0) + foci])[:, None]
    data = rng.poisson(data + 20).astype(np.uint16)
    info = ImageInfo("m.tif", "m", 1, 2, 1, 300, 300, "uint16", PX, None, "µm")
    return data, info, truth


def settings(**foci):
    s = config.Settings()
    s.segmentation.exclude_on_edges, s.segmentation.method = False, "threshold"
    s.segmentation.max_diameter = 100
    s.foci = [config.FociChannelSettings(channel=2, **foci)]
    return s


def test_scale_sizes():
    assert an.scale_sizes(0.5, 0.5) == [0.5]
    sizes = an.scale_sizes(0.5, 2.5)
    assert sizes[0] == pytest.approx(0.5) and sizes[-1] == pytest.approx(2.5)
    assert all(1.2 < b / a < 1.4 for a, b in itertools.pairwise(sizes))
    assert len(an.scale_sizes(0.05, 10)) == 8


def test_each_focus_once_at_its_own_size():
    data, info, truth = mixed()
    r = pipeline.analyze_frame(StageCache(), "m", data, info, settings(method="multi-scale", size=0.3, size_max=2.5))
    f = r.foci
    assert len(f) == len(truth)
    for sy, sx, size in truth:
        d = np.hypot(f.x / PX - sx, f.y / PX - sy)
        i = int(np.argmin(d))
        assert d.iloc[i] < 2.0
        assert size / 1.6 < f.detected_size.iloc[i] < size * 1.6  # found at about its own size


@pytest.mark.parametrize("dimension", ["3D", "max projection"])
@pytest.mark.parametrize("scope", ["per image", "per cell", "fixed value"])
def test_same_foci_as_watershed_on_even_foci(dimension, scope):
    data, info, spots = synthetic()
    s = config.Settings()
    s.segmentation.exclude_on_edges, s.segmentation.method = False, "threshold"
    s.foci = [config.FociChannelSettings(channel=2, method="multi-scale", size=0.5, size_max=1.5, threshold_scope=scope)]
    s.foci_general.dimension = dimension
    r = pipeline.analyze_frame(StageCache(), "k", data, info, s)
    assert sorted(r.cells["foci_ch2"]) == [3, 6]
    for _, sy, sx in spots:
        assert np.min(np.hypot(r.foci.x / PX0 - sx, r.foci.y / PX0 - sy)) < 1.0
