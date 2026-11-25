"""The update keeps the version that works until the new one is built."""
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("install", Path(__file__).parents[1] / "installer" / "install.py")
install = importlib.util.module_from_spec(spec)
spec.loader.exec_module(install)


def _old(target: Path) -> None:
    for name in ("program", ".venv"):
        (target / name).mkdir(parents=True)
        (target / name / "old.txt").write_text("old")


def test_a_failed_update_puts_the_old_version_back(tmp_path):
    _old(tmp_path)
    aside = install.set_aside(tmp_path, ["program", ".venv"])
    assert aside == ["program", ".venv"] and not (tmp_path / "program").exists()
    (tmp_path / "program").mkdir()  # the new version, half built
    (tmp_path / "program" / "new.txt").write_text("new")
    install.put_back(tmp_path, aside)
    for name in ("program", ".venv"):
        assert (tmp_path / name / "old.txt").is_file() and not (tmp_path / f"{name}.old").exists()
    assert not (tmp_path / "program" / "new.txt").exists()


def test_set_aside_is_all_or_none(tmp_path, monkeypatch):
    _old(tmp_path)
    real = Path.rename

    def rename(self, other):
        if self.name == ".venv":
            raise PermissionError("in use")
        return real(self, other)

    monkeypatch.setattr(Path, "rename", rename)
    with pytest.raises(SystemExit):
        install.set_aside(tmp_path, ["program", ".venv"])
    assert (tmp_path / "program" / "old.txt").is_file() and not (tmp_path / "program.old").exists()
