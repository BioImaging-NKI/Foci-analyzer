# Changes

Every entry says whether results change.

## 0.1.0 (in development)

First working version. There are no earlier results to compare with.

- Viewer: WebGL, tiles per pyramid level, smooth zoom, slice or max projection, contrast and colour per
  channel, outlines of cells and foci, markers that follow the threshold slider.
- Nuclei and cells: Cellpose, InstanSeg, StarDist, threshold, or labels from file; 2D or 3D.
- Foci: seeded watershed on a DoG, LoG or wavelet filtered image, or Spotiflow; extent at half maximum;
  threshold in units of the noise, per image, per cell or fixed.
- Colocalization per pair of foci channels: overlap or centroid distance, Pearson, Manders, chance control.
- Experiments with conditions, batch run, tables as CSV, Results page with tests (Holm corrected).
- Time-lapse: tracking of cells by overlap or Trackastra.
- Reads TIFF, CZI and LIF.
- Faster: the filters and the foci channels run at the same time on the processor; a batch with networks
  analyzes three images at a time; the networks are started when the window opens. Results do not change,
  except threshold segmentation (seeds found on a coarser grid; cell counts were the same on the test images).
- StarDist works, in an environment of its own (`stardist-env`, TensorFlow 2.15).
- Analyze page: nothing is analyzed until `Analyze this image` is pressed; `Live` then repeats it after every
  change. Cancel button in the status bar. Panel on the right with the channels (one slider with two ends,
  values, keys), the overlays per foci channel, the numbers of the image and the selected cell. Foci are
  lightly filled; colocalizing foci get a white ring and their share is shown. The projection the foci are
  detected on (max, EDF, sum) is shown in the viewer. No flicker when changing slices.
- Settings page: the networks, their versions, where they run (graphics card or processor) and from which
  Python; `Python of the models` moved there and can differ per network.
- Colocalization: a map in the viewer (white outline = colocalizing focus, white fill = overlap; key O),
  the share per pair of channels next to the image, and a Colocalization graph on the Results page with the
  chance control in grey. The chance control is now on by default (10 repeats) and in both directions: the
  cell table gets `chance_...` columns. `Colour cells by` fills the cells with a colour for any of their numbers.
- Large images (from 6 megapixels): `Visible part only` analyzes just what the viewer shows; a dashed frame
  marks that part. `Analyze all images` always uses whole images.
- StarDist is the default segmentation (it found the nuclei in all test sets and is the fastest). StarDist
  and InstanSeg say so when asked for 3D or (StarDist) a cell channel. Cellpose 3D on stacks of fewer than
  20 slices now segments every slice and links the slices (8 times faster, no over-splitting): results of
  3D Cellpose change.
- Wavelet filter: its noise estimate was too high, so it found far fewer foci than DoG at the same threshold;
  it now finds a comparable number. Results with the wavelet filter change.
- Spotiflow works on finely sampled images (it failed when the image had to be rescaled).
- Results page: `Tracks` graph (one line per tracked cell) and `Save graph data` (CSV of what the graph shows).
- `Also save label images`: a batch run can write the cells and the foci of every channel as label TIFFs.
- StarDist and InstanSeg in 3D ("2.5D"): every slice is segmented and the slices are linked where objects
  overlap (at least half of the smaller one); pieces of a nucleus that the network split in one slice stay one
  nucleus. On `3-channel` 24 nuclei (2D: 21); on the organoid 55 (2D: 26, nuclei on top of each other merge
  in a projection).
- Installers: `FociAnalyzer-Setup.exe` for Windows (choose the networks, or use your own environments) and
  `install.sh` for macOS and Linux. On Apple Silicon the PyTorch networks use the graphics chip.
  `installer/check.py` checks an installation. Results do not change.
- Foci no longer extend outside their cell (or outside the distance outside cells); in 3D this removed foci
  pixels above and below the cells. **Results change a little**: foci at a cell edge are a bit smaller.
- Cells segmented in 3D: foci are detected in 3D too (projections are greyed out); before, a projection
  could put foci where no cell is in that slice.
- A foci size smaller than a pixel no longer stops the analysis (IndexError).
- Smooth outlines, cell numbers, colour cells by label, an adjustable range for the colour scale, the
  selected cell highlighted, the distance outside cells as a grey line.
- Cells can be removed from the analysis by hand (key X); remembered in `removed_cells.tsv` in the results
  folder as points, so they stay out after a new analysis.
- Models chosen from a list per network (Cellpose 3: nuclei, cyto3, tissuenet_cp3; Cellpose 4: cpsam_v2,
  cpdino-vitb; StarDist; InstanSeg; Spotiflow) or a model of your own. Cellpose 4 has its own environment.
- Threshold segmentation: a slider for the threshold (× Otsu).
- Settings page: install a network with one click; use the graphics card or not; images at a time.
- Experiment page: stays open after opening images, updates the image list by itself, ✓ for analyzed images
  (also in the image list of the Analyze page). The status bar keeps its progress bar in one place.
- StarDist: when its nuclei come out smaller than 16 or larger than 120 pixels, the image is tried at four
  other scales and the result with nuclei closest to 40 pixels is kept (an uncalibrated high-resolution
  image gave 532 fragments, now 3 nuclei); *Expected diameter* now also enlarges. Results of the test
  images do not change. The Data Science Bowl 2018 model is no longer offered.
- Images for Fiji are saved as `.zip` (one uncompressed TIFF inside; Fiji opens it with the overlay), about
  a quarter of the size.
- Fixed: the dashed line at the distance outside cells was drawn around one group of cells only.
- Images analyzed in the viewer get a blue ✓; the progress bar sits left in the status bar.
- Fixed: an own Cellpose 4 environment given for Cellpose was reported as not installed, and its models
  greyed out. The automatic Cellpose model works with Cellpose 3 and 4.
- Windows: `Foci Analyzer.exe` in the installation folder, with the logo as icon (also for the shortcuts);
  Linux: the icon is also put in the icon theme.
- **Files named** matches with or without the extension (`*.tif` found nothing before).
- A user manual (`manual/`, not online yet).
- Logo in the window, the browser tab, the setup and the shortcuts.
- Colocalization map: a white line with a black edge, and overlap hatched black and white, visible on any colour.
- After `Analyze all images` the viewer shows the stored result of an image at once (same settings), also
  when it is opened from the Results page; contrast and channels shown stay when another image is opened.
- `Load settings…` and `Save settings…` (a results folder's settings.toml can be loaded too).
- Tables are tab-separated (`cells.tsv`, `foci.tsv`); older results folders with CSV are still read. `Save
  graph data` writes into the results folder.
- Save: the view (PNG), the whole image as shown (RGB TIFF), and the image for Fiji (TIFF with the outlines
  as overlay; also for all images with `Also save images with overlays`). In a tracked time-lapse the label
  images number a cell by its track.
- Status bar at the bottom: what the program is doing, step and progress bar, and when it is ready.
  Shortcut keys are shown next to their buttons; `[` and `]` step through slices when there is no time-lapse.
  Results do not change.
- Windows setup: the size of every network and the total for the choice; the page for your own
  environments only lists the networks left unticked (and is skipped when all are ticked); the installation's
  messages, step and progress bar are shown in the setup window (no console window); the Foci Analyzer logo
  in the wizard. Results do not change.
- Starting: the browser opens at once with a starting page, which turns into the window when the program is
  ready. On Windows a second Foci Analyzer (for example an installed one next to a development one) could
  take the same port, so the two answered at random; now the second start shows the first window.
- Setup: a Foci Analyzer still running from the installation folder is closed first (it made a reinstall
  fail); a failure shows the last messages, also when Python itself could not start
  (`install-output.log`). Cellpose is called Cellpose 3 there; the CellposeSAM download is only mentioned
  with Cellpose 4.
- Your own network environments may be conda environments: their library folders are put on PATH as
  `conda activate` would. Results do not change.
- Setup: a Python or environment left half-removed (an uninstall while the program ran) is removed and made
  again instead of failing; the folder choice is always shown.
- Updates: Foci Analyzer looks on GitHub for a newer version (at start, at most once a day; Settings page:
  on/off and Check for updates) and offers it with its release notes: Update now, Later, Not for a week,
  Skip this version. Updating downloads the new version, installs it with the same networks and starts it
  again; the browser tab follows. Results do not change.
- Setup: finds an existing installation (version, folder, networks) and offers Update/add networks or
  Remove; installed networks are shown as installed and only new ones count in the size. The window stays
  movable and scrollable while it installs (it was frozen).
- Cell border shift (was *Distance outside cells*): negative values move the border inwards, so foci near the
  edge of a cell are left out; the new checkbox *Foci in the whole image* replaces the old -1 (settings files
  with -1 are converted). **Results change** only when a negative shift is chosen.
- Several cells can be removed at once: Ctrl + click and Shift + drag a box select them, Esc clears.
- Viewer: the Z and T sliders are under the image, with a play button for time-lapses; Ctrl + wheel steps
  through Z, Alt + wheel through time. An image opened while the Analyze page was hidden no longer stays black.
- Overlays *Foci outlines* and *Foci markers* (were Foci and Markers); *Save images* with View, Full image and
  For Fiji, a green tick and Show in folder; *Colour cells by* in groups; the Threshold help explains the noise.
- The progress bar fills up when a task ends and then empties at once (it slid back); it no longer jumps back
  during an analysis. The starting page no longer reloads itself all the time (the logo blinked).
- Settings page: installing a network no longer asks the other networks for their version while it runs (red
  errors) and hides its messages when it succeeded.
- Setup: networks taken from your own environment are remembered and shown as such.
- Fixed: before the first analysis (or with Live off) only the first change of a setting was kept; later
  changes went into an old copy and were lost. Removed cells can be selected with a Shift box again (to put
  them back). With Live on, another time frame of an image analyzed in the viewer is analyzed too.
- *Cells to check*: a cell is flagged for saturated foci only when more than a tenth of its foci are saturated.
  *Colour cells by*: clearer names, and the channel in every foci entry.
