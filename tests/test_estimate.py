"""The foci size and threshold estimate on a made-up image with known foci."""

import numpy as np

from foci_analyzer import analysis as an


def test_estimate_on_synthetic_foci():
    rng = np.random.default_rng(3)
    img = rng.normal(100, 5, (1, 256, 256)).astype(np.float32)
    yy, xx = np.mgrid[:256, :256]
    sigma = 1.5  # px; FWHM 3.5 px = 0.7 µm at 0.2 µm pixels
    for _ in range(60):
        y, x = rng.uniform(20, 236, 2)
        img[0] += rng.uniform(150, 400) * np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / (2 * sigma**2))
    mask = np.ones((256, 256), bool)
    size = an.estimate_size(img[0], mask, 0.2)
    assert 0.5 < size < 1.0
    cal = an.Calibration(0.2, 0.2, "µm")
    s = an.sigma_px(size, cal, 1, 2.0)
    F = an.spot_filter(img, s, "DoG")
    unit, _ = an.noise_unit(F, np.ones((1, 256, 256), np.int32), "per image")
    _, h, _ = an.find_candidates(F / unit, s, 1.0, np.ones((1, 256, 256), np.int32), True)
    t, how = an.estimate_threshold(h)
    assert how in ("plateau", "otsu") and 40 <= int((h >= t).sum()) <= 75
