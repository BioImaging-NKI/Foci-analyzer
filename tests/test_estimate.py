"""The foci size and threshold estimate on made-up images with known foci."""

import numpy as np
import pytest

from foci_analyzer import analysis as an
from foci_analyzer.config import FociChannelSettings, FociGeneralSettings


def foci_image(n, fwhm_px, size=256, seed=3):
    rng = np.random.default_rng(seed)
    img = rng.normal(100, 5, (1, size, size)).astype(np.float32)
    s = fwhm_px / an.FWHM
    r = int(4 * s) + 1
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    for _ in range(n):
        y, x = rng.uniform(20, size - 20, 2)
        iy, ix = int(y), int(x)
        img[0, iy - r:iy + r + 1, ix - r:ix + r + 1] += rng.uniform(150, 400) * np.exp(
            -((yy - (y - iy)) ** 2 + (xx - (x - ix)) ** 2) / (2 * s**2))
    return img


def estimate(img, cells=None):
    cal = an.Calibration(0.2, 0.2, "µm")
    cellmap = np.ones(img.shape, np.int32) if cells is None else cells
    return an.estimate_watershed(img, cellmap, cal, FociGeneralSettings(), FociChannelSettings())


@pytest.mark.parametrize("n, size", [(20, 256), (60, 256), (2500, 1024)])
def test_estimate_whatever_the_number_of_foci(n, size):
    # FWHM 3.5 px = 0.7 µm at 0.2 µm pixels; few foci (mostly noise peaks) and very many (a tiled image)
    r = estimate(foci_image(n, 3.5, size))
    assert 0.55 <= r["size"] <= 0.85
    assert r["threshold"] > 2 and 0.7 * n <= r["foci"] <= 1.25 * n


def test_estimate_large_foci():
    r = estimate(foci_image(60, 7.0))  # 1.4 µm
    assert 1.2 <= r["size"] <= 1.6


def test_spread_sample_caps_per_cell():
    cells = np.r_[np.full(5000, 1), np.arange(2, 302)]  # one cell full of foci, 300 cells with one each
    pick = an.spread_sample(cells, 500, 10)
    assert len(pick) <= 500 and (cells[pick] == 1).sum() <= 10
    assert len(an.spread_sample(np.ones(2000, int), 500)) == 500  # one cell: still 500
