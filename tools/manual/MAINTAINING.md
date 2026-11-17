# Maintaining the user manual of Foci Analyzer

Notes for whoever (person or Claude, with the `update-manual` skill) updates the manual after a change to the
program. The manual: `manual/` (MkDocs Material, `mkdocs.yml` in the project root). **Not published yet**: it
goes on GitHub Pages once the repository has moved to BioImaging-NKI (workflow `manual.yml`, run by hand).
Rule: every noticeable change to the window is checked against the manual in the same piece of work.

## Tools

Documentation tools live in an environment of their own, next to the others (not in the project lock):
`C:/PythonProjects/Foci-Analyzer/docs-env` with mkdocs 1.6.1, mkdocs-material 9.8.1, mkdocs-glightbox 0.5.3,
playwright 1.65.2, pillow. Make it on another computer with
`uv venv docs-env --python 3.12` and `uv pip install --python docs-env/Scripts/python.exe mkdocs==1.6.1
mkdocs-material==9.8.1 mkdocs-glightbox==0.5.3 playwright==1.65.2 pillow`. Playwright drives the installed
Microsoft Edge (`channel="msedge"`): no browser download.

## Which code goes with which pages

| Changed code | Check these pages |
|---|---|
| `static/app.html`, `app.js`: Experiment page | `guide/experiment.md`, `start/first-analysis.md` (step 1) |
| `app.html`, `app.js`: settings panel, run box, Load/Save settings | `guide/analyze.md`, `start/first-analysis.md` (steps 3-5) |
| `app.html`, `app.js`, `viewer.js`: viewer, channels, overlays, selected cell, removal, Save | `guide/viewer.md`, `reference/keys.md` |
| keys (`keydown` in `app.js`) | `reference/keys.md`, `guide/viewer.md` |
| `static/results.js` | `guide/results.md` |
| Settings page (`loadProgram` in `app.js`, `server.program_info`) | `guide/settings-page.md` |
| `config.py` (labels, help, defaults) | regenerate `reference/settings.md`; `grep -rn "<label>" manual/` |
| `models.py` CATALOG / DOCS | regenerate `reference/models.md` |
| `pipeline.py`, `analysis.py` (methods) | `reference/analysis.md`, `help/faq.md` |
| table columns (`pipeline._tables`, `_cell_base`, `analysis.measure_foci`, `tracking.py`) | `reference/tables.md` |
| `batch.py`, `export.py`, `curation.py` (files written) | `guide/output-folder.md` |
| `installer/` | `start/install.md` |
| labels that README.md quotes | `README.md` |

Generated, never edit by hand: `manual/reference/settings.md` and `models.md` (`tests/test_manual.py` fails when
they are stale). Wrong text there: fix it in `config.py` or `models.py`.

## The example

- Folder `C:/PythonProjects/Foci-Analyzer/data/manual` (outside the synced folder): `images/foci_1.tif` and
  `foci_2.tif` (copies of `TestData/3-channel/Foci_02.tif` and `Foci_03.tif`: 1024 × 1024 × 8, 3 channels,
  0.24 µm pixels), `example.toml` (defaults with StarDist, nuclei channel 1, foci channels 2 and 3),
  `results/` (made by the screenshots), `home/` (the screenshot window's state; remade every time).
- The screenshot window makes one condition *Irradiated* with that folder (pattern `*.tif`) and results
  folder `results`, through the API.
- On another computer: copy the two images from TestData and remake `example.toml` (see `git log` of this file
  for the settings) — or move them into a release asset once the repository is public.

**Numbers quoted in the text** (recount after a change to the example or the analysis):

| Number | Page |
|---|---|
| foci_1 in the viewer: 21 nuclei, 1675 foci in channel 2, 1028 in channel 3 | `start/first-analysis.md` (step 4) |
| after Analyze all images: 21 and 26 nuclei, 79.8 and 90.2 foci of channel 2 per nucleus | `start/first-analysis.md` (step 5) |
| cell 17 selected, cell 18 removed in `viewer_cells.png` | `start/first-analysis.md` (caption) |

Recount with:

```python
import pandas as pd
c = pd.read_csv("C:/PythonProjects/Foci-Analyzer/data/manual/results/cells.tsv", sep="\t")
print(c.groupby("image").size().to_dict(), c.groupby("image")[["foci_ch2", "foci_ch3"]].mean().round(1).to_dict())
```

## Which screenshot shows what

| Screenshot (`manual/images/`) | Shows |
|---|---|
| `experiment.png` | the Experiment page with the example condition |
| `analyze.png` | the whole Analyze page after **Analyze this image** on foci_1 |
| `analyze_empty.png` | the viewer with foci_1 before analysis (not used in a page yet) |
| `settings_panel.png`, `display_panel.png` | the left and right panels of the Analyze page |
| `viewer_cells.png` | zoomed in on the largest nucleus (selected, yellow), its nearest neighbour removed (red), Numbers and Smooth outlines on |
| `selected_cell.png` | the Selected cell section of that nucleus |
| `results.png` | the Results page after **Analyze all images** |
| `settings_page.png` | the Settings page (lists this computer's graphics card and paths) |
| `logo.png` | copy of `src/foci_analyzer/static/logo.png` |

The removed cell is put back before **Analyze all images**, so `removed_cells.tsv` of the example is empty.

## Refresh and check

```powershell
C:\PythonProjects\Foci-Analyzer\docs-env\Scripts\python.exe tools/manual/refresh.py
```

Regenerates the reference pages, retakes all screenshots (own Foci Analyzer on port 8772, own home folder,
about 2 minutes with StarDist), builds the site strictly and runs `tests/test_manual.py`. Options:
`--shots analyze results` for some screenshots only; `--skip screenshots` for a change in the text only.

Then look at the changed images (`git status manual/images`) and browse the pages:
`docs-env/Scripts/python.exe -m mkdocs serve`.

## Publish (with the user's OK, after the move)

1. Move the repository to BioImaging-NKI; set `repo_url` (and `site_url`) in `mkdocs.yml`.
2. Settings → Pages → Source: GitHub Actions. Run the workflow `manual` by hand.
3. Check the live page itself (`curl -s <url> | grep -c "<img"`), not only that it loads.
