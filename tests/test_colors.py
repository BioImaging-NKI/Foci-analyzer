"""Channel colors: read from the file, fixed by channel number otherwise, foci in a lighter tint."""
import numpy as np
import tifffile

from foci_analyzer import export, io


def _lut(rgb):
    ramp = np.arange(256, dtype=np.uint8)
    return np.stack([ramp * (v // 255) for v in rgb]).astype(np.uint8)


def test_imagej_luts_are_read_and_gray_counts_as_none(tmp_path):
    path = tmp_path / "three.tif"
    data = np.zeros((3, 8, 8), np.uint16)
    tifffile.imwrite(path, data, imagej=True, metadata={"axes": "CYX", "mode": "composite",
                                                        "LUTs": [_lut((0, 0, 255)), _lut((255, 255, 255)), _lut((255, 0, 255))]})
    assert io.read_image(path)[1].colors == ["#0000ff", None, "#ff00ff"]
    blue = tmp_path / "blue.tif"  # every channel the same color: Foci Analyzer's colors instead
    tifffile.imwrite(blue, data, imagej=True, metadata={"axes": "CYX", "mode": "composite", "LUTs": [_lut((0, 0, 255))] * 3})
    assert io.read_image(blue)[1].colors is None
    plain = tmp_path / "plain.tif"
    tifffile.imwrite(plain, data)
    assert io.read_image(plain)[1].colors is None


def test_default_colors_do_not_depend_on_the_foci_channels():
    fa = export.default_colours(4, 1, 0)
    assert fa == [export.NUCLEI, export.OTHER[0], export.OTHER[1], export.OTHER[2]]
    assert export.default_colours(3, 1, 0, ["#ff0000", None, "#00ff00"]) == [(255, 0, 0), export.OTHER[0], (0, 255, 0)]


def test_foci_are_a_lighter_tint():
    assert export.foci_color((0, 255, 0)) == (128, 255, 128)
    assert export.foci_color((0, 0, 255)) == (128, 128, 255)
