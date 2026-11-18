"""Cell border shift: growing, shrinking (also between touching cells), the whole image, and old settings."""

import numpy as np

from foci_analyzer import analysis as an
from foci_analyzer import config


def _cal():
    return an.Calibration(pixel_size=1.0, z_spacing=1.0, unit="px")


def test_shrink_grow_and_whole_image():
    lab = np.zeros((1, 40, 60), np.uint16)
    lab[0, 10:30, 10:30] = 1
    lab[0, 10:30, 30:50] = 2  # touches cell 1
    area = lambda x, i: int((x == i).sum())
    small = an.expand(lab, _cal(), -3)
    assert area(small, 1) < area(lab, 1) and area(small, 2) < area(lab, 2)
    # shrunk away from the neighbor too: a gap between the two cells
    assert not small[0, 20, 27:33].any()
    big = an.expand(lab, _cal(), 3)
    assert area(big, 1) > area(lab, 1)
    everywhere = an.expand(lab, _cal(), np.inf)
    assert (everywhere > 0).all()


def test_old_minus_one_becomes_whole_image():
    s = config.from_dict({"foci_general": {"outside_distance": -1}})
    assert s.foci_general.whole_image and s.foci_general.outside_distance == 0
    s = config.from_dict({"foci_general": {"outside_distance": -1.5, "whole_image": False}})
    assert not s.foci_general.whole_image and s.foci_general.outside_distance == -1.5
    assert an.reach(config.from_dict({"foci_general": {"whole_image": True}}).foci_general) == np.inf
