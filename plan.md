# Foci Analyzer (Python) — plan

Python successor of the Fiji macro Foci Analyzer v1.96, with a browser GUI for a single local user.
It does not have to reproduce the Fiji numbers; it should be the better tool.

Status legend: **[done]** in the scaffold · **[next]** easy, no approval needed · **[plan]** difficult, wait for approval.

---

## 1. Goals and scope

In scope
- Nuclei/cell segmentation in 2D (on a projection) or 3D, independent of the foci dimension.
- Foci detection with (a) seeded watershed and (b) Spotiflow, for any number of foci channels.
- Foci detection in 3D, on a projection (max, EDF, sum) or per slice.
- Colocalization between foci channels, with more than one method and a chance control.
- Experiments: images from several folders, grouped in conditions (and replicates).
- A viewer with a tile pyramid, smooth zoom, and a slice / max-projection switch.
- Near real-time preview of foci detection with sliders; preview on the current view for large images.
- A results page with linked, clickable graphs. Everything stays viewable in the tool; table export is optional.
- Time-lapse as an optional part, with optional nucleus tracking.

Out of scope
- QuPath round trip, AreaMaxima detection, ImageJ `.zip` overlay images, bug-for-bug compatibility.

Bugs of v1.96 that disappear by design: StarDist rescale precedence, Z centroid offset, 2D labels used on 3D
statistics, foci split over two nuclei, colocalization on undilated nuclei, broken "load settings" keys,
3D size filter ignoring Z spacing, infinite Z kernels for 2D images.

---

## 2. Technical choices

| Topic | Choice | Why |
|---|---|---|
| Python | 3.12, `uv`, environment outside surfdrive (`C:/PythonProjects/Foci-Analyzer/.venv`) | Same as FAST-HIPPOS; best support by torch-based tools |
| Server | Standard-library `ThreadingHTTPServer` on 127.0.0.1, JSON + binary endpoints, server-sent events for progress | Single user, no build step, no extra dependencies, same pattern as FAST-HIPPOS |
| Frontend | Plain JavaScript modules, no Node, no bundler | Same reason; easy to maintain |
| Viewer | WebGL2, raw 16-bit tiles as textures, compositing in a shader | Contrast, channel toggles and zoom cost nothing; see §6 |
| Graphs | Plotly (one bundled file, no Node), as in Microirradiation-recruitment-analysis; WebGL traces when there are many cells | Click, lasso, zoom and SVG/PNG export are built in, so far less code to write than own canvas plots |
| Core numerics | numpy, scipy, scikit-image on the CPU | Fast enough for preview (see §5.6); GPU only for deep-learning models |
| Deep learning | PyTorch only: Cellpose, Spotiflow, Trackastra, (InstanSeg) | One framework; StarDist needs TensorFlow, so it is an optional extra |
| Settings | Dataclasses with metadata → schema → form | One definition gives the GUI, validation, TOML files and the manual |
| Storage | A project folder: `project.toml`, per image a folder with labels (TIFF), tables (Parquet) and settings | Results reopen in the tool; export to TSV/XLSX on request |

Decision needed (D1): deep-learning models in this environment (extras) or in your existing environments
`C:/PythonProjects/Cellpose3`, `Cellpose4`, `Spotiflow` through a subprocess, as FAST-HIPPOS does for Cellpose.
Proposal: support both; extras by default, an "other Python" path under Advanced.

---

## 3. Data model

```
Experiment
 ├─ conditions[]          name, colour, folders[], file pattern, replicate rule
 ├─ images[]              path, series, condition, replicate, status, per-image overrides
 ├─ settings              one Settings object for the whole experiment (§4)
 └─ results/              per image: labels, foci tables, QC; experiment tables
```

- A condition has one or more folders and a file pattern. Replicate comes from a sub-folder name, a file-name
  pattern, or is set by hand. **[done]** conditions and scanning; **[next]** replicates, multi-series files.
- Tables (one row per object, tidy):
  - `cells`: image, condition, replicate, frame, cell id, track id, size, shape, per-channel intensities,
    per-foci-channel count and summed measures, colocalization measures.
  - `foci`: image, cell id, channel, position (µm and px), size, peak, mean, integrated intensity above local
    background, signal-to-noise, probability (Spotiflow), distance to the cell edge, colocalization partner.
- A focus belongs to exactly one cell: the cell under its centroid (or the nearest cell within the allowed
  distance outside). No splitting.
- Internal array order is always `T, C, Z, Y, X`; calibration in µm travels with the array. **[done]**

Readers: TIFF/OME-TIFF with `tifffile` **[done]**; `.czi`, `.lif`, `.nd2` through `bioio` plug-ins as extras
**[next]**; large files read lazily per plane **[plan, with the viewer]**.

---

## 4. Settings and GUI layout

One browser tab with three pages:

1. **Experiment** — conditions, folders, the image list with status. **[done: basic]**
2. **Analyze** — settings on the left, viewer in the middle, live counts on the right. **[done: shell]**
3. **Results** — statistics (§8). **[plan]**

Settings sections; each shows 2–5 settings, the rest sits under a collapsed *Advanced*. **[done]**

| Section | Visible | Advanced |
|---|---|---|
| Image | nuclei channel, cell channel | XY binning, crop border, axis order |
| Segmentation | method, 2D/3D, min and max diameter (µm) | model, diameter, thresholds, projection used, exclude edges, labels folder |
| Foci (general) | detect in 3D / max / EDF / sum / per slice | distance outside cells |
| Foci channel (repeatable) | channel, method, foci size (µm), sensitivity | min/max size, extent rule, threshold scope, fixed threshold, Spotiflow model |
| Colocalization | on/off, method | min overlap, max distance, chance control |
| Time-lapse | tracking method | linking distance, gap closing |

Foci size is in µm (pixels when the image has no calibration). The form shows the value in pixels next to it.

---

## 5. Analysis design

The analysis is a chain of cached stages. Each stage has a key made of its own settings and the key of its
input, so a changed slider only recomputes what depends on it. **[done: cache; stages follow]**

```
read → prepare (bin, crop, projection) → segment → filter foci channel → candidates
     → accept (threshold) → extents → assign to cells → measure → colocalize
```

### 5.1 Preparation **[next]**
Binning, crop, projections (max, sum; EDF by local variance focus measure with a smoothed height map).
Segmentation and foci detection choose their dimension separately; a 2D label image is extended over Z when
foci are detected in 3D. That combination stays the default for stacks with few slices.

### 5.2 Segmentation **[plan]**
One interface, `segment(image, settings) → labels`, with these back-ends:

| Back-end | Use | Remark |
|---|---|---|
| Cellpose 3 (`nuclei`, `cyto3`) / 4 (`cpsam`) | nuclei or cells, 2D, 2.5D, 3D | `cpsam` is slow on a 6 GB GPU; Cellpose 3 stays supported |
| StarDist `2D_versatile_fluo` | nuclei 2D | optional extra (TensorFlow); rescale to 0.25 µm/px as in Fiji |
| InstanSeg (to evaluate) | nuclei 2D | PyTorch, light; candidate default if it holds up on your data |
| Threshold + watershed | nuclei 2D | no model needed; always available |
| Labels from file | any | label image or ImageJ ROI zip |

After any back-end: remove on edges (XY), size filter with the real voxel size, relabel. Manual removal by
clicking cells in the viewer, stored with the image.
Decision needed (D2): the default back-end. Proposal: decide after a comparison on your data (Cellpose 3
`nuclei`, StarDist, InstanSeg) in milestone 3.

### 5.3 Foci detection: seeded watershed, revised **[plan]**

Points where I would deviate from the macro, with the reason:

1. **Filter in physical units.** A difference of Gaussians with σ from the foci size in µm; in Z the σ follows
   from the Z spacing and a fixed axial elongation (default 3×), not from an empirical square root.
2. **Background without the outlier filter.** The wide Gaussian is taken from a foci-suppressed image (grey
   opening with a radius of about one foci diameter, on a down-sampled copy for speed). Same intent as
   "remove outliers", fewer halos around bright foci, no dependence on integer rounding.
3. **Threshold in noise units of the filtered image.** The macro measures the SD in the raw image and applies it
   to the filtered image, so the meaning of "3×" changes with the foci size. New: noise σ = 1.4826 × MAD of the
   filtered values inside the cells; the slider is a signal-to-noise level (log scale). MAD barely reacts to the
   foci themselves, so cells full of foci do not raise their own threshold.
4. **Threshold scope is explicit**: per image (default), per cell, or one fixed value for the whole experiment.
   A fixed value is the honest choice when conditions are compared and the control has no foci; the GUI says so.
5. **Detection and extent are separate.** Candidates are local maxima of the filtered image; each has a peak
   height. *Accepting* a candidate is a comparison with the threshold. The *extent* of a focus is its watershed
   basin, cut either at the threshold (as in Fiji) or at half of its own peak (size no longer depends on the
   threshold). Decision needed (D3): default extent rule. Proposal: half-maximum.
6. **Measurements above local background**: integrated intensity = Σ(raw − local background) over the extent,
   next to mean, peak and signal-to-noise. "Total foci intensity" per cell becomes a real integral.

### 5.4 Foci detection: Spotiflow **[plan]**
- Run once per image at a low probability; keep position, probability and intensity per spot.
- The slider filters on probability (and optionally intensity) without running the model again.
- Spots are points. For size and overlap an extent is derived as in §5.3 point 5, seeded by the spots (optional).
- 2D model on projections or per slice; the 3D model for stacks when available for the data type.
- Models: `general` by default, own models by path.

### 5.5 Colocalization **[plan]**
For every pair of foci channels (not only A and B):
- **Object overlap**: a focus colocalizes when its overlap with a focus of the other channel is ≥ n voxels or
  ≥ a fraction of its own size. Per focus: partner id and overlap fraction.
- **Centroid distance**: nearest neighbour in the other channel within d µm (anisotropy included); works for
  Spotiflow points.
- **Per cell**: count and fraction colocalized for both directions, Manders M1/M2 on the foci masks
  (intensity-weighted), Pearson inside the cell as an extra.
- **Chance control**: expected colocalization when channel B foci are placed at random positions inside the same
  cell (n repeats). Reported next to the observed value, so a dense channel is not mistaken for colocalization.

### 5.6 Real-time preview **[plan]**
- Candidates carry their peak height (or probability), so the threshold slider is applied **in the browser**:
  markers and counts per cell change while dragging, without a request.
- The server recomputes extents after the slider rests (debounced), for the current view only when the image
  is large. This replaces the Alt-ROI trick: "preview on the current view" is a toggle.
- Changing the foci size recomputes the filter stage: about 0.1–0.3 s for 1024 × 1024 × 8 on the CPU (to be
  measured in milestone 2). Segmentation is never recomputed by foci sliders.
- "Apply to all images" runs the batch in a background thread with progress and cancel.

### 5.7 Time-lapse and tracking **[plan, last]**
- Frames are a dimension, not split files. Segmentation and foci run per frame.
- Tracking links existing labels: built-in overlap/distance assignment (no extra dependency) and Trackastra as
  an optional extra (handles divisions). Ultrack solves segmentation and tracking together and needs more
  set-up; it is left out unless Trackastra falls short.
- Results: track id per cell, traces per cell on the Results page, lineage for divisions.
- All time-lapse settings and views stay hidden when the images have one frame.

---

## 6. Viewer **[plan]**

Server
- `GET /api/image/region?id&t&c&z|mip&level&x&y&w&h` returns raw 16-bit pixels. **[done]**
- Pyramid levels by 2× averaging, made on first request and cached. **[done, in memory]**
- Label tiles (cells, foci) as 32-bit; outlines are drawn in the shader from label changes, so they stay one
  screen pixel wide at any zoom.
- Large images: planes read lazily from the file and levels cached on disk.

Browser
- WebGL2: one texture per channel and tile, a fragment shader applies min/max, colour and additive mixing.
- Zoom glides to its target (exponential approach, as in FAST-HIPPOS); the coarser level stays on screen until
  the finer tiles have arrived.
- Slice / max-projection switch; in projection mode the foci of all slices are shown, in slice mode the foci
  near the slice (marker size shows the Z distance).
- Layers: channels, cell outlines, foci extents (semi-transparent), foci markers (drawn from the table, so
  the slider can filter them directly), cell numbers and counts as text.
- Click a cell: its row in the tables and its point in the graphs. Click a point in a graph: the viewer goes
  to that cell.

The scaffold has a simple canvas viewer (whole plane, one level) that already uses the region endpoint, the
slice / projection switch and smooth wheel zoom. The WebGL viewer replaces it.

---

## 7. Batch run and outputs **[plan]**

- Run all or the selected images; resume after a stop; images with results are skipped unless settings changed.
- Per image a folder: `labels_cells.tif`, `labels_foci_ch*.tif`, `cells.parquet`, `foci.parquet`,
  `settings.toml`, `qc.json`. Experiment level: combined tables.
- Export on request: TSV or XLSX, label images, viewer snapshots, graphs as SVG/PNG.

---

## 8. Results page **[plan]**

- **Per-cell distributions by condition**: strip/violin/box, with replicate means on top ("superplot").
- **Tests on replicate means**, not on cells, to avoid inflated significance; the page states n per level.
- **Scatter** of any two columns (density when crowded), **histogram/ECDF**, **foci per cell** bar/ECDF,
  **colocalization** observed versus chance.
- **Image QC table**: cells, foci per cell, threshold, noise, flags (saturation, too few cells, drift of
  threshold between images).
- **Time-lapse**: traces per cell, mean ± band per condition.
- Tests can use cells or experiment means as replicates (the user chooses, as in
  Microirradiation-recruitment-analysis); p values are Holm-corrected; small groups are reported as such.
- **Cells to check**: quality flags (for example saturated, on the size limit, outlier in its condition).
  Flags never exclude a cell; the user decides.
- Results that are older than the current settings are marked as out of date.
- Selection is shared: lasso or click in one graph highlights the same cells elsewhere; double click opens the
  cell in the viewer. Filters (for example cell size, intensity gate for transfected cells) apply to all graphs.

---

## 9. Milestones

| # | Content | State |
|---|---|---|
| 0 | Scaffold: package, settings schema, TIFF reader, experiment scan, server, GUI shell, simple viewer, tests | done |
| 1 | WebGL tile viewer with label and marker layers | done |
| 2 | Foci by seeded watershed (§5.3) with live preview (§5.6); threshold segmentation and labels from file | done |
| 3 | Deep-learning segmentation back-ends | done: Cellpose, InstanSeg, StarDist (tested; StarDist in its own `stardist-env`) |
| 4 | Spotiflow | done |
| 5 | Measurements, colocalization, batch run, storage | done |
| 6 | Results page | done (first version) |
| 7 | Time-lapse and tracking | done: overlap linking and Trackastra, tested on a made-up time-lapse only |
| 8 | Manual (MkDocs), example data, desktop shortcut | installer and shortcut done; manual not started |

### As built: where it differs from the plan above

- **Storage**: `cells.csv`, `foci.csv`, `settings.toml`, `run.json` in the results folder (no Parquet: one
  dependency less, and the tables open in Excel). Label images are not saved; a click on a result opens the
  image and repeats the (cached) analysis.
- **Images in memory**: every opened image is read completely (largest test image: 7533 × 7506 × 4 channels,
  450 MB). The viewer does request tiles per pyramid level, so the browser stays light. Reading regions
  straight from pyramidal CZI files is not built.
- **Live preview**: the threshold slider filters the candidates in the browser at once (markers and count);
  the outlines follow when the server has redone the extents (0.1 s in 2D, about 0.7 s for 1024² × 8 in 3D).
  "Current view only" for large images is not built; a whole 22-megapixel 2D image takes 1.3 s per step.
- **Networks**: they run in a worker process with its own Python (`models-env` next to the app environment,
  or any Python given in the settings), so an existing Cellpose environment can be used. Images sampled
  finer than 0.25 µm per pixel are made smaller for the networks.
- **Spotiflow on stacks** with a 2D model: spots come from the max projection, their slice from the
  brightest value. A 3D model (name with a 3D model) is used directly.
- **Filters**: DoG (default), LoG and à trous wavelet; all 3D for 3D detection, sigma in µm per axis with an
  axial elongation factor (default 2).
- **Default threshold** 8 × noise (robust SD inside the cells). On `TestData/3-channel` 10 × gives about the
  Fiji count in 3D; projections need about 5 ×. It stays a visual choice, hence the slider.
- **Batch**: images are analyzed in parallel processes when no network is used (8 test images in 20 s);
  with networks one after another, because they share the graphics card.
- **Extra**: Quit button and "start twice opens the running window", needed without a console window.

### Measured on this computer (Quadro RTX 3000)

| Step | 1024² × 8, 2 foci channels | 5679 × 3852, 1 channel |
|---|---|---|
| Cellpose (cyto3/nuclei, automatic diameter) | 5.6 s | about 45 s |
| InstanSeg | 4.1 s | about 10 s |
| Filter + candidates + extents (CPU) | 4 s | 5 s |
| One slider step | 0.7 s | 1.3 s |
| Model loading, first use per session | 25–30 s per network | |

InstanSeg found the nuclei of the cultured cells in `3-channel` (28, Cellpose 24 of which edge nuclei are
removed) and is clearly faster on large images; Cellpose was more robust on the image without a real
nuclear stain (`Hi-res 3D`). Default: Cellpose.

### Open, in order of use to you

1. **Speed**. Done: filters and foci channels in parallel threads (foci part of a 1024² × 8 image with two
   channels 2.9 → 1.9 s), batch with networks three images at a time, networks started when the window
   opens, faster threshold segmentation (1.7 → 0.7 s); the whole 3D test image 4.5 → 2.8 s, as max
   projection 2.1 → 1.1 s. Measured and not done: filters on the GPU (the
   processor part is now 1 to 2 s per image, and sending the image to the model process costs about as much
   as it would win). "Visible part only" for images from 6 megapixels: the preview analyzes what the viewer
   shows (at least 1024 px wide), a dashed frame marks it.
2. **Results page**: done: colocalization graph with chance control, tracks graph, `Save graph data`.
   Open: more flags.
3. **Manual** with `/anthropic-skills:update-manual` (methods are already in the context help).
4. Cellpose 4 (`cpsam`) through an existing environment (Settings page > own Python for cellpose):
   written, not tested.
5. **Foci tracking over time** (registration of the nucleus, then linking foci): foci have ids and
   positions per frame and a `track` column of their cell, nothing more yet.
6. Comparison with the Fiji outputs (foci per cell, StarDist nuclei): `3-channel` Fiji 41 and 26 (ch 2, 3),
   ours 43 and 37 at threshold 12 in 3D; `Deconvolved with cytoplasm` Fiji 10.5, ours 4 at foci size 0.5 µm
   and 8 to 12 at 0.2 µm; `Hi-res 3D` Fiji 259 (ch 1), ours 120 at 0.5 µm and 250 at 0.3 µm. The foci size
   setting decides the count as much as the threshold does: no defect found, not tuned further.
7. Tested on all test sets with every method (2026-11-16): StarDist found nuclei everywhere and is the
   default; Cellpose failed on `1-channel` (0 cells) and `Hi-res 3D` first image (1 cell); threshold fails on
   the dense organoid (1 cell). Not tested: Cellpose 4, the installer's StarDist step, real time-lapse data.
8. Not built: foci tracking, the manual. Cellpose 4 (`cpsam`) was tested through the existing environment
   `C:/PythonProjects/Cellpose4/.venv` (works; about 28 s per 1024² image on the Quadro RTX 3000).

---

## 10. Way of working

- Tokens matter more than building speed: one session, no `/orchestrate`, sub-agents only when they clearly
  save reading.
- `CHANGELOG.md`: every entry says whether results change.
- Result files are saved atomically (write to a temporary file, then rename).
- Commits per milestone, pushed to the private repository `bvandenbroek/Foci-Analyzer-Python` (to be moved
  to BioImaging-NKI).

## 11. Decisions

- **D1** Networks in a separate environment by default; the installer asks, and any existing environment
  can be given in the settings. Done.
- **D2** Default segmentation: Cellpose (see the measurements above); InstanSeg for speed.
- **D3** Foci extent at half maximum. Done.
- **D4** Git repository with a commit per milestone. Done.
- **D5** Package and command: `foci_analyzer` / `foci-analyzer`.
