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


def _grid(flip_x: bool, swap: bool):
    """A 3 x 2 Leica tile scan of 100 x 80 px tiles at 1 µm per pixel, 10% overlap, stage axes as a microscope may have
    them: X mirrored and/or X and Y swapped."""
    tiles = []
    for fy in range(2):
        for fx in range(3):
            x, y = fx * 90e-6 + 0.01, fy * 72e-6 + 0.02
            x = -x if flip_x else x
            x, y = (y, x) if swap else (x, y)
            tiles.append({"FieldX": str(fx), "FieldY": str(fy), "PosX": str(x), "PosY": str(y)})
    return tiles


@pytest.mark.parametrize(("flip_x", "swap"), [(False, False), (True, False), (False, True), (True, True)])
def test_lif_tiles_go_where_the_stage_was(flip_x, swap):
    off = io.tile_offsets(_grid(flip_x, swap), 1.0, 100, 80)
    assert off.tolist() == [[0, 0], [90, 0], [180, 0], [0, 72], [90, 72], [180, 72]]


def test_lif_tiles_without_positions_go_on_the_grid():
    tiles = [{"FieldX": str(i % 2), "FieldY": str(i // 2)} for i in range(4)]
    assert io.tile_offsets(tiles, 1.0, 100, 80).tolist() == [[0, 0], [100, 0], [0, 80], [100, 80]]


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
