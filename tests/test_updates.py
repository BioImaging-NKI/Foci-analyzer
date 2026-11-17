"""The update check: versions, and which release counts (only those with the setup attached)."""

import json

from foci_analyzer import updates


def test_version_order():
    k = updates.version_key
    assert k("0.1.0.dev0") < k("0.1.0") < k("v0.1.1") < k("0.2.0") < k("1.0") < k("2.0.0rc1") < k("2.0.0")


def test_latest_skips_fiji_drafts_and_prereleases(tmp_path, monkeypatch):
    setup = [{"name": "FociAnalyzer-Setup-0.3.0.exe", "browser_download_url": "https://x/s.exe"}]
    releases = [
        {"tag_name": "v1.97", "assets": [{"name": "Foci_Analyzer.ijm"}], "body": "Fiji"},
        {"tag_name": "v0.2.0", "assets": setup, "body": "two", "zipball_url": "z2", "html_url": "h2"},
        {"tag_name": "v0.4.0", "assets": setup, "draft": True},
        {"tag_name": "v0.5.0", "assets": setup, "prerelease": True},
        {"tag_name": "v0.3.0", "assets": setup, "body": "three", "zipball_url": "z3", "html_url": "h3",
         "published_at": "2026-12-01T10:00:00Z"},
    ]
    f = tmp_path / "releases.json"
    f.write_text(json.dumps(releases), encoding="utf-8")
    monkeypatch.setenv("FOCI_ANALYZER_UPDATE_URL", f.as_uri())
    new = updates.latest()
    assert new["version"] == "0.3.0" and new["notes"] == "three" and new["source"] == "z3" and new["date"] == "2026-12-01"


def test_preferences_snooze_and_skip(tmp_path, monkeypatch):
    monkeypatch.setattr(updates, "__version__", "0.1.0")
    p = updates.Preferences(tmp_path / "updates.json")
    p.data["latest"] = {"version": "0.2.0"}
    assert p.show()
    p.data["skip"] = "0.2.0"
    assert not p.show()
    p.data.update(skip="", snooze_until=4e12)
    assert not p.show()
    p.save()
    assert updates.Preferences(tmp_path / "updates.json").data["snooze_until"] == 4e12
