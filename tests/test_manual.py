"""The manual stays in line with the program: the generated reference pages are up to date, every page in the
navigation exists, and every image a page shows is there."""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUAL = ROOT / "manual"
sys.path.insert(0, str(ROOT / "tools" / "manual"))


def test_reference_pages_are_current():
    import gen_reference
    assert (MANUAL / "reference" / "settings.md").read_text(encoding="utf-8") == gen_reference.render(), \
        "run: python tools/manual/gen_reference.py"
    assert (MANUAL / "reference" / "models.md").read_text(encoding="utf-8") == gen_reference.render_models(), \
        "run: python tools/manual/gen_reference.py"


def test_pages_and_images_exist():
    nav = re.findall(r":\s+([\w/-]+\.md)\s*$", (ROOT / "mkdocs.yml").read_text(encoding="utf-8"), re.MULTILINE)
    assert nav, "no pages in mkdocs.yml"
    for page in nav:
        path = MANUAL / page
        assert path.is_file(), page
        for img in re.findall(r"\]\(([^)\s]+\.png)\)", path.read_text(encoding="utf-8")):
            assert (path.parent / img).resolve().is_file(), f"{page}: {img}"
