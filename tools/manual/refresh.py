"""Bring the manual in line with the program: regenerate the reference pages, retake the screenshots, build the
site strictly and run the manual tests (tools/manual/MAINTAINING.md).

    <docs-env python> tools/manual/refresh.py [--shots analyze results] [--skip screenshots]
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
APP_PYTHON = os.environ.get("FA_PYTHON", "C:/PythonProjects/Foci-Analyzer/.venv/Scripts/python.exe")


def step(title: str, cmd: list[str]) -> None:
    print(f"\n== {title}", flush=True)
    r = subprocess.run(cmd, cwd=ROOT, check=False)
    if r.returncode:
        raise SystemExit(f"failed: {title}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--shots", nargs="*", default=[], help="only the screenshots whose names start with these")
    p.add_argument("--skip", default="", help="comma-separated: screenshots, build, tests")
    a = p.parse_args()
    skip = set(a.skip.split(",")) - {""}
    step("reference pages", [APP_PYTHON, "tools/manual/gen_reference.py"])
    if "screenshots" not in skip:
        step("screenshots", [sys.executable, "tools/manual/screenshots.py", *a.shots])
    if "build" not in skip:
        os.environ["NO_MKDOCS_2_WARNING"] = "1"
        step("strict build", [sys.executable, "-m", "mkdocs", "build", "--strict", "-q"])
    if "tests" not in skip:
        step("manual tests", [APP_PYTHON, "-m", "pytest", "-q", "-p", "no:warnings", "tests/test_manual.py"])
    print("\nDone. Look at the changed screenshots: git status manual/images")


if __name__ == "__main__":
    main()
