# Foci Analyzer (Python) — notes for Claude Code

Python successor of the Fiji macro Foci Analyzer (v1.96, `C:\Fiji\plugins\Foci-Analyzer\Foci_Analyzer.ijm`,
same as GitHub BioImaging-NKI/Foci-analyzer). Browser window for a single local user. It does not have to
reproduce the Fiji results. Design, state per milestone and open items: `plan.md` (keep it current, §9).

## Environment
- The project folder is in surfdrive (synced). The uv environment is OUTSIDE it: always
  `export UV_PROJECT_ENVIRONMENT=C:/PythonProjects/Foci-Analyzer/.venv` before `uv` commands, or call
  `C:/PythonProjects/Foci-Analyzer/.venv/Scripts/python.exe` directly. The file `.venv` in the project folder
  is a placeholder that stops uv from making an environment in surfdrive: never delete it.
- Networks: `C:/PythonProjects/Foci-Analyzer/models-env` (torch 2.6 cu124, cellpose 3.1, spotiflow,
  instanseg-torch, trackastra), `cellpose4-env` (cellpose 4.2.1.1; Cellpose 3 and 4 cannot share one) and
  `stardist-env` (Python 3.11, TensorFlow 2.15, stardist; TensorFlow 2.21 in models-env was 100 times slower,
  so keep them apart). `models.CATALOG` lists the pretrained models per network and which environment runs
  them; `models.installed()` reads the package folders (else asks that Python); keys `cellpose` (any version: the
  automatic model), `cellpose3`, `cellpose4`. A user's Cellpose 4 given for `cellpose` also serves `cellpose4`. Found automatically because they are next to `.venv`.
  Per computer another Python can be given on the Settings page (`~/.foci_analyzer/program.json`).
- Tests: `python -m pytest -q -p no:warnings`; lint: `python -m ruff check src tests`.
- The user's window state is in `~/.foci_analyzer`. Test with a separate one:
  `.claude/launch.json` → `foci-analyzer-test` (port 8771, `--home C:/PythonProjects/Foci-Analyzer/test-home`).
  Python changes need a restart of that server; files in `static/` are read on every request.
  In the test window `window.FA = { viewer, V, S }` gives the state.
- No Node, no conda. GPU: Quadro RTX 3000, 6 GB. The user's own environments: `C:/PythonProjects/Cellpose3`,
  `Cellpose4`, `Spotiflow`.
- Git: private repository `BioImaging-NKI/Foci-analyzer` (moved there on 2026-11-18; the Fiji macro is now
  `BioImaging-NKI/Foci-analyzer-Fiji`; BioImaging-NKI is a user account, bvandenbroek has push rights); commit and push
  per milestone is agreed. `TestData/` and `LiFT/` stay out of it.
- After every change in `static/*.js`: load the test window and check that `window.FA` exists before
  committing. A `//` comment put in the middle of a one-line statement swallows the rest of the line and
  the whole window stays empty (this happened twice); use `/* */` there.
- The project path contains `&`: in cmd it must be quoted; write multi-line Python to a file in the
  scratchpad instead of a heredoc, and edit project files with the Edit/Write tools, not `sed` with `\n`.

## Modules (`src/foci_analyzer/`)
- `config.py`: settings as dataclasses made with `opt()`; `schema()` drives the form, `from_dict()` validates.
  A new setting is added there only (label, context help, `advanced`, `slider`, `when`). At most 5 visible
  settings per section.
- `io.py`: TIFF, CZI (scenes), LIF (series) → array `T, C, Z, Y, X` + `ImageInfo`. `experiment.py`:
  conditions (folders or files) → image list; a file with several series gives several images (`id`).
- `analysis.py`: the image analysis as plain functions; all arrays `Z, Y, X` (Z = 1 for 2D).
- `pipeline.py`: one frame as cached stages: prepare → segment → per foci channel filter → candidates →
  extents → cells → tables. `stages.py`: the cache, with time per stage.
- `models.py` + `model_worker.py`: networks in a worker process with its own Python (one per environment).
  The worker must not import `foci_analyzer`. `App.cancel()` kills the workers; `warm_up()` starts them early.
- Preview: `App.frame()` computes and stores the result in `App.last`; overlays and candidates only read
  `App.last` (they never start an analysis). `StageCache.now` feeds the status bar, `.cancelled` stops it.
- `export.py`: the store the viewer reads back after a batch run (`.viewer/*.npz` in the results folder:
  labels and candidates per frame), label TIFFs, RGB rendering, ImageJ TIFFs with outline overlays (roifile).
  Overlay TIFFs must stay uncompressed (they are zipped instead: Fiji opens a .zip with one TIFF, overlay
  intact) and their ROIs must have no `group` and no label options: Fiji drops
  the overlay of a compressed TIFF, colours ROIs by group, and draws numbers otherwise (checked with Fiji).
  Fiji can be asked from here: `C:/Fiji/fiji-windows-x64.exe --headless --run script.py` (Jython; write
  results to a file, the console stays empty; Flatten does not work headless).
- `curation.py`: cells removed by hand, as points in `removed_cells.tsv` in the results folder; applied to
  every result (viewer and batch) so the tables leave them out.
- `batch.py`: all images (parallel processes without networks), tables to the results folder.
  `tracking.py`: overlap linking, or Trackastra. `stats.py`: tests with Holm correction.
- `cli.py` takes the port and opens the browser on a starting page before the slow imports, then hands the
  server to the app (`LocalServer`: no SO_REUSEADDR on Windows, else two programs share a port).
- `updates.py`: the update check (GitHub releases of `BioImaging-NKI/Foci-analyzer` that have the setup
  attached: the Fiji releases live there too) and the update (source archive -> `installer/install.py
  --update`, which waits for the program, serves the progress on its port and starts the new version).
  `FOCI_ANALYZER_UPDATE_URL` points the check at a test JSON (file:// works). Choices in `~/.foci_analyzer/updates.json`.
- `server.py`: standard-library server on 127.0.0.1; `/api/*` needs the token written into the page and a
  local Host header; only images of the experiment can be read.
- `static/`: `app.js` (pages, form, preview, batch), `viewer.js` (WebGL2 tile viewer), `results.js`
  (Plotly graphs, tests), `util.js`; plain JavaScript modules, no build step. Plotly comes from the
  `plotly` Python package.
- `installer/`: `install.py` (the installation on every system, standard library + uv; `--networks-only`
  is what the Install buttons of the Settings page run), `branded_launcher` copies uv's
  `foci-analyzer-window.exe` to `Foci Analyzer.exe` with the logo (the data uv appends after the program part
  is kept), `FociAnalyzer.iss`
  + `build-windows.ps1` (setup.exe; Inno Setup 6 is installed per user here; the wizard shows the progress
  that `install.py --progress FILE` writes), `run-install.cmd`, `check.py`
  (every installed network on a made-up image). `install.sh` (macOS, Linux) and `install.bat` (developers)
  call `install.py`. Network versions are pinned in `install.py`; change them only after testing them together.
- `.github/workflows/build.yml`: tests on three systems per push; installers only by hand or for a tag (the
  private repository pays for minutes, macOS counts ten times).

## Working agreements
- Present a detailed plan before implementing difficult parts.
- Keep the number of visible settings low; sections with a collapsed Advanced part; context help everywhere.
- Tokens matter more than building speed: no `/orchestrate`, few sub-agents, read only what is needed.
- Analysis speed matters: parallel where possible, GPU where it helps.
- The manual: `manual/` (MkDocs Material, `mkdocs.yml`), kept with `/anthropic-skills:update-manual`; notes in
  `tools/manual/MAINTAINING.md` (page map, example, screenshots, refresh). Not online until the move to
  BioImaging-NKI. Tools in `C:/PythonProjects/Foci-Analyzer/docs-env`; `tests/test_manual.py` checks the
  generated reference. A change to the window means checking the manual in the same piece of work.
- Inspiration (same author, both readable with `gh`; Microirradiation is private):
  Jalink-lab/FAST-HIPPOS-Python and Jalink-lab/Microirradiation-recruitment-analysis.
- gh CLI: `"C:/Program Files/GitHub CLI/gh.exe"` (not on PATH in Git Bash).
- Organization default: metric units. American spelling (color, gray, center, canceled) in code, window and manual.
- Test data: `TestData/` (channels as far as known): `3-channel` and `1-channel`
  (nuclei 1, foci 2 and 3), `Hi-res 3D` (no nuclear stain; nuclei from channel 2), `Deconvolved with
  cytoplasm` (float32, no Z spacing), `Multi-series` (LIF, 5 series), `Organoid` (CZI stack),
  `Spheroid pyramidal` (large 2D CZI). Fiji outputs are in sub-folders; they are good, not ground truth.
