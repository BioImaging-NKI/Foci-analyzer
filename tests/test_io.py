"""Reading images: calibration of OME-TIFF files, color (RGB) CZI files."""
import numpy as np
import pytest
import tifffile

from foci_analyzer import io


def test_ome_tiff_calibration(tmp_path):
    path = tmp_path / "ome.ome.tif"
    tifffile.imwrite(path, np.zeros((4, 8, 8), np.uint16), ome=True,
                     metadata={"axes": "ZYX", "PhysicalSizeX": 65.0, "PhysicalSizeXUnit": "nm", "PhysicalSizeY": 65.0,
                               "PhysicalSizeYUnit": "nm", "PhysicalSizeZ": 0.3, "PhysicalSizeZUnit": "µm"})
    info = io.read_image(path)[1]
    assert info.unit == "µm" and info.pixel_size == pytest.approx(0.065) and info.z_spacing == pytest.approx(0.3)


def test_rgb_czi_gives_three_channels(tmp_path):
    czi = pytest.importorskip("pylibCZIrw.czi")
    path = tmp_path / "rgb.czi"
    img = np.zeros((20, 30, 3), np.uint8)
    img[..., 0], img[..., 1], img[..., 2] = 10, 20, 30  # stored as blue, green, red
    with czi.create_czi(str(path)) as w:
        w.write(data=img, plane={"C": 0, "T": 0, "Z": 0})
    arr, info = io.read_image(path)
    assert arr.shape == (1, 3, 1, 20, 30) and list(arr[0, :, 0, 0, 0]) == [30, 20, 10]
    assert info.colors == ["#ff0000", "#00ff00", "#0000ff"]
