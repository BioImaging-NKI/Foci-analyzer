# Foci Analyzer (Python) — notes for Claude Code

Python successor of the Fiji macro Foci Analyzer (v1.96, `C:\Fiji\plugins\Foci-Analyzer\Foci_Analyzer.ijm`,
same as GitHub BioImaging-NKI/Foci-analyzer). Browser window for a single local user. It does not have to
reproduce the Fiji results. Design, milestones and open decisions: `plan.md` (keep it current).

## Environment
- The project folder is in surfdrive (synced). The uv environment is OUTSIDE it: always
  `export UV_PROJECT_ENVIRONMENT=C:/PythonProjects/Foci-Analyzer/.venv` before `uv` commands, or call
  `C:/PythonProjects/Foci-Analyzer/.venv/Scripts/python.exe` directly. The file `.venv` in the project folder
  is a placeholder that stops uv from making an environment in surfdrive: never delete it.
- Tests: `python -m pytest -q -p no:warnings`; lint: `python -m ruff check src tests`.
- The user's window state is in `~/.foci_analyzer`. Test with a separate one:
  `.claude/launch.json` → `foci-analyzer-test` (port 8771, `--home C:/PythonProjects/Foci-Analyzer/test-home`,
  conditions point at `TestData/`).
- No Node, no conda. GPU: Quadro RTX 3000, 6 GB. Existing environments that may be reused for models:
  `C:/PythonProjects/Cellpose3`, `Cellpose4`, `Spotiflow`.
- Not a git repository yet (plan.md, D4). Commit and push only when the user says so.

## How it works now (milestone 0, scaffold)
- `config.py`: settings as dataclasses made with `opt()`; `schema()` drives the form, `from_dict()` validates.
  A new setting is added there only. At most 5 visible settings per section; the rest `advanced=True`.
- `io.py`: TIFF → array `T, C, Z, Y, X` + `ImageInfo` (µm, or `px` when not calibrated).
- `experiment.py`: conditions (name, folders, pattern) → image list.
- `imagestore.py`: opened images, planes (slice or max projection) at pyramid levels, in memory.
- `stages.py`: `StageCache`, the basis for the cached analysis chain (plan.md §5). No analysis stages yet.
- `server.py`: standard-library server on 127.0.0.1; `/api/*` needs the token written into the page and a
  local Host header; only images of the experiment can be read.
- `static/app.{html,css,js}`: pages Experiment | Analyze | Results; plain JavaScript modules, no build step.
  The viewer is a simple 2D-canvas stand-in for the WebGL tile viewer of plan.md §6.

## Working agreements
- Present a detailed plan before implementing difficult parts; milestones 1 and 2 wait for approval.
- Keep the number of visible settings low; sections with a collapsed Advanced part.
- Tokens matter more than building speed: no `/orchestrate`, few sub-agents, read only what is needed.
- Inspiration (same author, both readable with `gh`; Microirradiation is private):
  Jalink-lab/FAST-HIPPOS-Python (standard-library server, plain JavaScript, tile pyramid viewer, TOML settings)
  and Jalink-lab/Microirradiation-recruitment-analysis (Setup | Review | Results, conditions, Plotly graphs,
  statistics per cell or on experiment means with Holm correction, quality flags that never exclude).
- gh CLI: `"C:/Program Files/GitHub CLI/gh.exe"` (not on PATH in Git Bash).
- Organization default: metric units.
- Test data: `TestData/` (1024², 8-bit, 3 channels × 8 slices, 0.240 µm XY, 0.755 µm Z; and a 1-channel
  version). More data and Fiji outputs will come from the user.
