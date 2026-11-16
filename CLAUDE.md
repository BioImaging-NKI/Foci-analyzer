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
  instanseg-torch, trackastra) and `stardist-env` (Python 3.11, TensorFlow 2.15, stardist; TensorFlow 2.21 in
  models-env was 100 times slower, so keep them apart). Found automatically because they are next to `.venv`.
  Per computer another Python can be given on the Settings page (`~/.foci_analyzer/program.json`).
- Tests: `python -m pytest -q -p no:warnings`; lint: `python -m ruff check src tests`.
- The user's window state is in `~/.foci_analyzer`. Test with a separate one:
  `.claude/launch.json` → `foci-analyzer-test` (port 8771, `--home C:/PythonProjects/Foci-Analyzer/test-home`).
  Python changes need a restart of that server; files in `static/` are read on every request.
  In the test window `window.FA = { viewer, V, S }` gives the state.
- No Node, no conda. GPU: Quadro RTX 3000, 6 GB. The user's own environments: `C:/PythonProjects/Cellpose3`,
  `Cellpose4`, `Spotiflow`.
- Git: private repository `bvandenbroek/Foci-Analyzer-Python` (to move to BioImaging-NKI); commit and push
  per milestone is agreed. `TestData/` and `LiFT/` stay out of it.
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
- `batch.py`: all images (parallel processes without networks), tables to the results folder.
  `tracking.py`: overlap linking, or Trackastra. `stats.py`: tests with Holm correction.
- `server.py`: standard-library server on 127.0.0.1; `/api/*` needs the token written into the page and a
  local Host header; only images of the experiment can be read.
- `static/`: `app.js` (pages, form, preview, batch), `viewer.js` (WebGL2 tile viewer), `results.js`
  (Plotly graphs, tests), `util.js`; plain JavaScript modules, no build step. Plotly comes from the
  `plotly` Python package.
- `install.bat`: installation for users without Python (uv, app environment, optional networks, shortcut).

## Working agreements
- Present a detailed plan before implementing difficult parts.
- Keep the number of visible settings low; sections with a collapsed Advanced part; context help everywhere.
- Tokens matter more than building speed: no `/orchestrate`, few sub-agents, read only what is needed.
- Analysis speed matters: parallel where possible, GPU where it helps.
- The manual is made with `/anthropic-skills:update-manual` (not started; not online yet).
- Inspiration (same author, both readable with `gh`; Microirradiation is private):
  Jalink-lab/FAST-HIPPOS-Python and Jalink-lab/Microirradiation-recruitment-analysis.
- gh CLI: `"C:/Program Files/GitHub CLI/gh.exe"` (not on PATH in Git Bash).
- Organization default: metric units.
- Test data: `TestData/` (channels as far as known): `3-channel` and `1-channel`
  (nuclei 1, foci 2 and 3), `Hi-res 3D` (no nuclear stain; nuclei from channel 2), `Deconvolved with
  cytoplasm` (float32, no Z spacing), `Multi-series` (LIF, 5 series), `Organoid` (CZI stack),
  `Spheroid pyramidal` (large 2D CZI). Fiji outputs are in sub-folders; they are good, not ground truth.
