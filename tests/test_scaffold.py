"""Tests of the scaffold: settings, reading, experiment scan, stages and the server."""

import json
import threading
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pytest
import tifffile

from foci_analyzer import config, experiment
from foci_analyzer.imagestore import downsample
from foci_analyzer.io import read_image
from foci_analyzer.stages import StageCache

TESTDATA = Path(__file__).parent.parent / "TestData"


# ---------------------------------------------------------------- settings
def test_settings_roundtrip(tmp_path):
    s = config.Settings()
    s.foci.append(config.FociChannelSettings(channel=3, method="spotiflow"))
    s.segmentation.min_diameter = 6.0
    config.save(s, tmp_path / "s.toml")
    assert config.load(tmp_path / "s.toml") == s


def test_settings_validation():
    with pytest.raises(ValueError, match="not one of"):
        config.from_dict({"segmentation": {"method": "magic"}})
    with pytest.raises(ValueError, match="smaller"):
        config.from_dict({"foci": [{"channel": 0}]})
    with pytest.raises(ValueError, match="at least one"):
        config.from_dict({"foci": []})
    assert config.from_dict({"image": {"unknown": 1, "nuclei_channel": 2}}).image.nuclei_channel == 2
    assert config.from_dict({"foci": [{"size": 1}]}).foci[0].size == 1.0  # a whole number for a float setting


def test_schema_has_every_setting_and_hides_most():
    for sec in config.schema():
        assert sec["fields"] and all(f["label"] for f in sec["fields"])
        methods = next((f["choices"] for f in sec["fields"] if f["key"] == "method"), [None])
        for m in methods:  # what is shown at the same time: at most five settings outside Advanced
            shown = [f for f in sec["fields"] if not f["advanced"] and (not f["when"] or m in f["when"].get("method", [m]))]
            assert len(shown) <= 5, (sec["key"], m)


# ---------------------------------------------------------------- reading
def _write(path, arr, axes, **meta):
    tifffile.imwrite(path, arr, imagej=True, resolution=(1 / 0.25, 1 / 0.25), metadata={"axes": axes, **meta})


def test_read_orders_axes_and_calibration(tmp_path):
    arr = np.arange(2 * 3 * 4 * 5, dtype=np.uint16).reshape(3, 2, 4, 5)  # Z, C, Y, X
    _write(tmp_path / "a.tif", arr, "ZCYX", unit="micron", spacing=0.7)
    data, info = read_image(tmp_path / "a.tif")
    assert data.shape == (1, 2, 3, 4, 5)
    assert np.array_equal(data[0, 1, 2], arr[2, 1])
    assert info.unit == "µm" and info.pixel_size == pytest.approx(0.25) and info.z_spacing == pytest.approx(0.7)


def test_read_without_calibration(tmp_path):
    tifffile.imwrite(tmp_path / "b.tif", np.zeros((6, 7), np.uint8))
    data, info = read_image(tmp_path / "b.tif")
    assert data.shape == (1, 1, 1, 6, 7) and info.unit == "px" and info.pixel_size is None


@pytest.mark.skipif(not (TESTDATA / "3-channel" / "Foci_02.tif").exists(), reason="no test data")
def test_read_test_data():
    data, info = read_image(TESTDATA / "3-channel" / "Foci_02.tif")
    assert data.shape == (1, 3, 8, 1024, 1024) and info.unit == "µm"
    assert info.pixel_size == pytest.approx(0.2403, abs=1e-3) and info.z_spacing == pytest.approx(0.7553, abs=1e-3)
    data1, _ = read_image(TESTDATA / "1-channel" / "Foci_02_channel2.tif")
    assert data1.shape == (1, 1, 8, 1024, 1024)  # 8 slices, not 8 channels
    assert np.array_equal(data1[0, 0], data[0, 1])


def test_downsample_keeps_edge_pixels():
    plane = np.arange(15, dtype=np.float32).reshape(3, 5)
    small = downsample(plane, 1)
    assert small.shape == (2, 3)
    assert small[0, 0] == pytest.approx(plane[:2, :2].mean()) and small[1, 2] == pytest.approx(plane[2, 4])


# ---------------------------------------------------------------- experiment
def test_scan_conditions(tmp_path):
    for folder, names in {"ctrl": ["a.tif", "b.tif", "notes.txt"], "treated": ["c.tiff"], "treated/sub": ["d.tif"]}.items():
        (tmp_path / folder).mkdir(parents=True, exist_ok=True)
        for n in names:
            (tmp_path / folder / n).write_bytes(b"")
    conds = experiment.conditions_from_dicts([
        {"name": "control", "folders": [str(tmp_path / "ctrl")]},
        {"name": "treated", "folders": [str(tmp_path / "treated"), str(tmp_path / "ctrl"), str(tmp_path / "gone")],
         "recursive": True}])
    images, remarks = experiment.scan(conds)
    assert [(i.name, i.condition) for i in images] == [("a", "control"), ("b", "control"), ("c", "treated"), ("d", "treated")]
    assert any("folder not found" in r for r in remarks) and sum("is in both" in r for r in remarks) == 2
    with pytest.raises(ValueError, match="two conditions"):
        experiment.conditions_from_dicts([{"name": "x"}, {"name": "x"}])


# ---------------------------------------------------------------- stages
def test_stage_cache_recomputes_only_downstream():
    cache = StageCache()

    def chain(sigma, threshold):
        k1, a = cache.run("filter", {"sigma": sigma}, ("image",), lambda: sigma * 10)
        return cache.run("accept", {"threshold": threshold}, (k1,), lambda: a + threshold)[1]

    assert chain(1, 5) == 15 and cache.computed == ["filter", "accept"]
    assert chain(1, 6) == 16 and cache.computed == ["filter", "accept", "accept"]  # the filter was reused
    assert chain(2, 6) == 26 and cache.computed[-2:] == ["filter", "accept"]  # both again


# ---------------------------------------------------------------- server
@pytest.fixture
def server(tmp_path, monkeypatch):
    monkeypatch.setenv("FOCI_ANALYZER_HOME", str(tmp_path / "home"))
    from foci_analyzer.server import make_server
    srv, app = make_server(0)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}", app
    srv.shutdown(); srv.server_close()


def _call(base, path, token, body=None):
    req = urllib.request.Request(base + path, data=json.dumps(body).encode() if body is not None else None,
                                 headers={"X-Token": token, "Content-Type": "application/json"})
    with urllib.request.urlopen(req) as r:
        return r.read(), r.headers


def test_server_end_to_end(server, tmp_path):
    base, app = server
    page = urllib.request.urlopen(base + "/").read().decode()
    assert app.token in page and "{{TOKEN}}" not in page
    with pytest.raises(urllib.error.HTTPError) as e:  # no token, no access
        _call(base, "/api/state", "wrong")
    assert e.value.code == 403

    state = json.loads(_call(base, "/api/state", app.token)[0])
    assert [s["key"] for s in state["schema"]][:2] == ["image", "segmentation"]

    with pytest.raises(urllib.error.HTTPError) as e:
        _call(base, "/api/settings", app.token, {"settings": {"segmentation": {"method": "magic"}}})
    assert e.value.code == 400 and "not one of" in json.loads(e.value.read())["error"]
    state["settings"]["foci"][0]["size"] = 0.8
    saved = json.loads(_call(base, "/api/settings", app.token, {"settings": state["settings"]})[0])
    assert saved["foci"][0]["size"] == 0.8
    assert config.load(tmp_path / "home" / "last_settings.toml").foci[0].size == 0.8

    arr = np.random.default_rng(0).integers(0, 4000, (2, 3, 9, 11), dtype=np.uint16)  # C, Z, Y, X
    (tmp_path / "data").mkdir()
    _write(tmp_path / "data" / "im.tif", arr.transpose(1, 0, 2, 3), "ZCYX", unit="micron", spacing=1.0)
    exp = json.loads(_call(base, "/api/experiment", app.token,
                           {"conditions": [{"name": "c", "folders": [str(tmp_path / "data")]}]})[0])
    path = exp["images"][0]["path"]
    q = urllib.parse.quote(path)

    info = json.loads(_call(base, f"/api/image/info?path={q}", app.token)[0])
    assert (info["channels"], info["slices"], info["width"]) == (2, 3, 11) and len(info["ranges"]) == 2
    body, headers = _call(base, f"/api/image/region?path={q}&c=1&z=mip", app.token)
    got = np.frombuffer(body, "<u2").reshape(int(headers["X-Height"]), int(headers["X-Width"]))
    assert np.array_equal(got, arr[1].max(axis=0))
    body, headers = _call(base, f"/api/image/region?path={q}&c=0&z=2&x=3&y=1&w=4&h=5", app.token)
    assert np.array_equal(np.frombuffer(body, "<u2").reshape(5, 4), arr[0, 2, 1:6, 3:7])
    body, headers = _call(base, f"/api/image/region?path={q}&c=0&z=0&level=1", app.token)
    assert (headers["X-Height"], headers["X-Width"]) == ("5", "6")

    for bad in (f"/api/image/region?path={q}&c=7&z=0", f"/api/image/info?path={urllib.parse.quote(str(tmp_path / 'x.tif'))}"):
        with pytest.raises(urllib.error.HTTPError) as e:
            _call(base, bad, app.token)
        assert e.value.code in (400, 404)
