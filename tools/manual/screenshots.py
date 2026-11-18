"""Screenshots of the user manual (manual/images), of the example (MAINTAINING.md, "The example"), taken with
the installed Microsoft Edge driven by Playwright. A Foci Analyzer of its own is started for this (own home
folder, port 8772), so the state of the user's window is not touched; no file dialogs are used.

    <docs-env python> tools/manual/screenshots.py [only the shots whose names start with these, e.g. analyze]

Environment: FA_PYTHON (the Python of Foci Analyzer), FA_EXAMPLE (the example folder).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

from PIL import Image
from playwright.sync_api import Page, sync_playwright

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "manual" / "images"
PYTHON = os.environ.get("FA_PYTHON", "C:/PythonProjects/Foci-Analyzer/.venv/Scripts/python.exe")
EXAMPLE = Path(os.environ.get("FA_EXAMPLE", "C:/PythonProjects/Foci-Analyzer/data/manual"))
HOME = EXAMPLE / "home"
PORT = 8772
URL = f"http://127.0.0.1:{PORT}/"
VIEW = {"width": 1440, "height": 900}


def save(target, name: str, **kw) -> None:
    path = OUT / f"{name}.png"
    target.screenshot(path=path, animations="disabled", **kw)
    im = Image.open(path)
    im.save(path, optimize=True)
    print(f"{name}.png {im.size[0]}×{im.size[1]}")


def home_crop(name: str, x0: int = 0) -> None:
    """A 16:10 part of a window screenshot, from x0 below the bar: for the cards of the home page."""
    im = Image.open(OUT / f"{name}.png")
    w = min(im.size[0] - x0, 1400)
    im.crop((x0, 80, x0 + w, 80 + w * 10 // 16)).save(OUT / f"home_{name}.png", optimize=True)
    print(f"home_{name}.png")


def call(page: Page, path: str, body: dict | None = None):
    """An API request from inside the page (it holds the token)."""
    return page.evaluate("""async ([path, body]) => {
      const r = await fetch(path, { method: body ? 'POST' : 'GET', body: body ? JSON.stringify(body) : undefined,
        headers: { 'X-Token': window.TOKEN, 'Content-Type': 'application/json' } });
      const j = await r.json(); if (!r.ok) throw new Error(j.error); return j; }""", [path, body])


def start() -> subprocess.Popen:
    if HOME.exists():
        shutil.rmtree(HOME)  # a fresh window: no remembered display or experiment
    results = EXAMPLE / "results"
    for name in ("removed_cells.tsv",):
        (results / name).unlink(missing_ok=True)
    proc = subprocess.Popen([PYTHON, "-m", "foci_analyzer", "--port", str(PORT), "--no-browser", "--home", str(HOME)],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for _ in range(100):
        try:
            urllib.request.urlopen(URL, timeout=1).read(100)
            return proc
        except OSError:
            time.sleep(0.2)
    proc.kill()
    raise SystemExit("Foci Analyzer did not start")


def setup(page: Page) -> None:
    """The example experiment and settings, as a user would have made them on the Experiment page."""
    call(page, "/api/settings/load", {"path": str(EXAMPLE / "example.toml")})
    call(page, "/api/experiment", {"conditions": [{"name": "Irradiated", "folders": [str(EXAMPLE / "images")],
                                                   "pattern": "*.tif", "recursive": False}],
                                   "output": str(EXAMPLE / "results")})


def wait_preview(page: Page) -> None:
    page.wait_for_function("window.FA && window.FA.V.preview && document.getElementById('status').textContent.startsWith('Ready')",
                           timeout=300_000)
    page.wait_for_timeout(1500)  # the tiles and outlines


def main(want_prefixes: list[str]) -> None:
    want = (lambda n: not want_prefixes or any(n.startswith(p) for p in want_prefixes))
    OUT.mkdir(parents=True, exist_ok=True)
    proc = start()
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(channel="msedge")
            page = browser.new_page(viewport=VIEW, device_scale_factor=1.5, color_scheme="light")
            page.goto(URL)
            page.wait_for_function("window.FA")
            setup(page)
            page.goto(URL)  # the page reads the settings and the experiment when it opens
            page.wait_for_function("window.FA && window.FA.S.images.length === 2")
            if want("experiment"):
                save(page, "experiment")
                home_crop("experiment", 240)  # the page is centered

            # Analyze: the first image, analyzed in the viewer
            page.click("[data-page=analyze]")
            page.select_option("#image-select", index=1)
            page.wait_for_timeout(2500)
            if want("analyze_empty"):
                save(page.locator("#viewer"), "analyze_empty")
            page.click("#analyze-now")
            wait_preview(page)
            if want("analyze"):
                save(page, "analyze")
            if want("settings_panel"):
                save(page.locator("#side"), "settings_panel")
            if want("display_panel"):
                save(page.locator("#display"), "display_panel")

            # zoomed in on one nucleus: numbers, smooth outlines, a selected and a removed cell
            page.check("#show-numbers")
            page.check("#smooth")
            page.wait_for_timeout(1500)
            # the largest nucleus is selected; its nearest neighbor is removed (and put back afterwards)
            sel, rem = page.evaluate("""() => { const ps = [...FA.viewer.numbers.polys], area = (p) => (p.bb[2] - p.bb[0]) * (p.bb[3] - p.bb[1]);
              const big = ps.sort((a, b) => area(b) - area(a))[0];
              const near = ps.filter((p) => p !== big).sort((a, b) => Math.hypot(a.cx - big.cx, a.cy - big.cy) - Math.hypot(b.cx - big.cx, b.cy - big.cy))[0];
              return [[big.cx, big.cy], [near.cx, near.cy]]; }""")
            page.evaluate("([x, y]) => FA.viewer.onclick(x, y)", rem)
            page.click("#toggle-cell")
            page.wait_for_timeout(1200)
            page.evaluate("([x, y]) => { FA.viewer.onclick(x, y); FA.viewer.centerOn(x, y, 300); }", sel)
            page.wait_for_timeout(2500)
            if want("viewer_cells"):
                save(page.locator("#canvas-wrap"), "viewer_cells")
            if want("selected_cell"):
                page.locator("#cell-info").scroll_into_view_if_needed()
                save(page.locator("#cell-info").locator("xpath=.."), "selected_cell")
            page.evaluate("([x, y]) => FA.viewer.onclick(x, y)", rem)
            page.click("#toggle-cell")  # back in the analysis, so that the batch run counts every cell
            page.wait_for_timeout(800)
            page.uncheck("#show-numbers")
            page.uncheck("#smooth")

            # Foci detection (guide/foci.md): a foci channel's settings, and a few nuclei raw, filtered and at a
            # low threshold, with the foci markers
            if want("foci_channel"):
                group = page.locator(".group", has=page.locator(".estimate")).first
                group.evaluate("(g) => g.scrollIntoView({ block: 'start' })")  # above the run panel at the bottom
                page.wait_for_timeout(300)
                save(group, "foci_channel")
            if any(want(n) for n in ("foci_raw", "foci_filtered", "foci_low")):
                # only foci channel 2 (channel 1 is the nuclei, channel 3 the other foci), on a nucleus near the
                # middle of the image, without a selected cell
                page.evaluate("() => FA.viewer.onclick(1, 1)")
                page.uncheck("#ch0")
                page.uncheck("#ch2")
                page.uncheck("#show-foci")  # the markers only
                page.uncheck("#show-coloc")
                page.locator("#foci-toggles input").nth(1).uncheck()
                page.evaluate("""() => { const W = FA.V.info.width / 2, H = FA.V.info.height / 2;
                  const p = [...FA.viewer.numbers.polys].sort((a, b) => Math.hypot(a.cx - W, a.cy - H) - Math.hypot(b.cx - W, b.cy - H))[0];
                  FA.viewer.centerOn(p.cx, p.cy, 220); }""")
                page.wait_for_timeout(2500)
                if want("foci_raw"):
                    save(page.locator("#canvas-wrap"), "foci_raw")
                page.click("#filter-mode [data-f=filtered]")
                page.wait_for_timeout(2500)
                if want("foci_filtered"):
                    save(page.locator("#canvas-wrap"), "foci_filtered")
                page.click("#filter-mode [data-f=raw]")
                threshold = page.locator('input[type=range][aria-label="Threshold"]').first
                threshold.evaluate("(r) => { r.value = 3; r.dispatchEvent(new Event('input', { bubbles: true })); }")
                page.wait_for_timeout(1500)
                wait_preview(page)
                if want("foci_low"):
                    save(page.locator("#canvas-wrap"), "foci_low")
                threshold.evaluate("(r) => { r.value = 8; r.dispatchEvent(new Event('input', { bubbles: true })); }")
                page.wait_for_timeout(1500)
                wait_preview(page)
                page.check("#ch0")
                page.check("#ch2")
                page.check("#show-foci")
                page.check("#show-coloc")
                page.locator("#foci-toggles input").nth(1).check()

            # Analyze all images, then the Results page
            page.click("#run")
            page.wait_for_function("document.querySelector('#page-results') && !document.querySelector('#page-results').hidden",
                                   timeout=900_000)
            page.wait_for_timeout(3000)
            if want("results"):
                save(page, "results")
                home_crop("results")
            if want("settings_page"):
                page.click("[data-page=settings]")
                page.wait_for_selector("#program table", timeout=120_000)
                page.wait_for_timeout(800)
                save(page, "settings_page")
            call(page, "/api/quit", {})
        print(json.dumps({"images": str(OUT)}))
    finally:
        time.sleep(0.5)
        proc.kill()


if __name__ == "__main__":
    main(sys.argv[1:])
