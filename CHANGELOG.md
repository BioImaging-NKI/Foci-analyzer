# Changes

Every entry says whether results change.

## Not released yet

- Updates: a failed update puts the version that was there back (the old program and its environment are set
  aside until the new ones are installed), and a release candidate now counts as newer than a beta of the
  same version (`2.0.0rc1` after `2.0.0b5`).
- OME-TIFF files without ImageJ metadata are calibrated: pixel size, Z spacing and frame interval from the
  OME data (they opened in pixels). Results change for such files: sizes in µm instead of pixels, and the
  µm settings apply as meant.
- Color (RGB) CZI files open as three channels, red, green and blue (only the blue part was read). Results
  change for such files.
- LIF tile scans that were not stitched on the microscope open as one image, every tile at its stage position
  (only the first tile was read). Results change for such files.
- An image without a pixel size gets a remark that the µm settings are read as pixels.

## 2.0.0b5 (2026-11-23)

- The update check in the background runs every 6 hours instead of once a day, also while Foci Analyzer
  stays open, so a new version shows up without **Check for updates**.
- The viewer shows its zoom (*Zoom 66%*), and **100%** next to **Fit** shows one image pixel on one screen
  pixel. Results do not change.
- Smooth outlines follow the nuclei and cells more closely (the spline cut the corners of elongated and
  dented nuclei). Results do not change.
- Channel colors: **Colors** above the channels takes the colors stored in the image file (ImageJ TIFF,
  OME-TIFF, CZI, LIF; not when they give every channel the same color, as some LIF files do), or Foci
  Analyzer's own. Foci Analyzer's colors now go by channel number, so changing the foci channels no longer
  changes the colors. Foci are drawn in a lighter tint of their channel color,
  which stands out on their own signal (also in the saved pictures and Fiji overlays). The colored bars of
  the foci channels in the settings follow a change of channel. Results do not change.
- Colocalization map: where foci of two channels overlap is filled in the mix of their colors (green and red
  give yellow) instead of black and white stripes, and only from a zoom of 150%. Results do not change.

## 2.0.0b4 (2026-11-22)

- **Segment nuclei** / **Segment cells** (key S, left of **Analyze this image**) finds only the nuclei or
  cells, to try the segmentation settings quickly. **Estimate** works on them, and **Analyze this image**
  afterwards does not segment again. **Smooth outlines** moved to key E.
- **Visible part only** sits right of **Analyze this image** and works for segmenting too. Nuclei and cells
  found in a visible part are kept: as long as the segmentation settings stay the same, a part inside an
  earlier one takes its cells from there instead of segmenting again.
- **Estimate** (watershed) measures the foci size on the foci above a threshold it finds itself, sampled
  over all cells, instead of on the 300 brightest peaks: it now works for images with few foci and for large
  tiled images with many. Results change only when **Estimate** is used.
- Smooth outlines of nuclei and cells are fitted with a spline (in the viewer and in the Fiji overlays); a
  nucleus or cell in several pieces now shows all of them (one piece could vanish). Foci are drawn in the
  color of their channel. Results do not change.
- **Full screen** at the top right (or F11) uses the whole screen; Esc or F11 goes back.
- Only the newest browser tab of Foci Analyzer works; older tabs say so and can take over again.
- The starting page: a thicker progress ring in the colors of the manual's title, the logo centered in it.

## 2.0.0b3 (2026-11-20)

- Windows setup: the installation stopped at *[1/4] The program* with *Failed to create Python minor version
  link directory ... untrusted mount point (os error 448)*. Inno Setup 6.7 passes Windows' RedirectionGuard
  on to the installation, which then cannot use the link uv makes in the Python folder; the setup now turns it
  off. Results do not change.

## 2.0.0b2 (2026-11-19)

- New foci method **multi-scale**, for foci of very different sizes: the filter at several sizes between the
  smallest and the largest foci size (a slider with two ends), each focus at the size where it stands out
  most, only round peaks, one focus per large focus. The size at which each focus was found is the new column
  `detected_size` of `foci.tsv`. **Estimate** proposes the size range too. Results of the other methods do not
  change.
- Watershed: **Round foci only** (Advanced, off by default) leaves out peaks that are not round, such as the
  edge of a bright nucleus or a side peak of a larger focus (1 to 8 % of the foci on test images). Results
  change only when it is turned on.
- The starting page shows how far the start is: a ring around the logo and the current step, timed as at the
  last start.
- Spotiflow works with any *Foci size* and with the 3D models: the image is now resized by Foci Analyzer
  itself (Spotiflow's own resizing refused 3D models and most sizes, with *3D scaling is not supported yet*
  or *Invalid scale factor*). Results change only where these errors occurred.
- Windows: `Foci Analyzer.exe` has the logo again (newer uv launchers keep their start data inside the program
  file, which the setup did not expect).
- **Estimate** has a **?** like the settings. Manual: a page *Foci detection* on the methods, every foci
  setting and how to tune them. Results do not change.

## 2.0.0b1 (2026-11-18)

The first public version of Foci Analyzer in Python, the successor of the Fiji macro (1.x, now in
[Foci-analyzer-Fiji](https://github.com/BioImaging-NKI/Foci-analyzer-Fiji)). A beta: please report problems.
It does not reproduce the Fiji results exactly; it has its own methods.

- A window in the browser on your own computer: experiments with conditions, a live preview of the analysis
  in a viewer for large 2D and 3D images and time-lapses, and a Results page with graphs and statistics.
- Nuclei and cells with StarDist, Cellpose 3 or 4, InstanSeg, a threshold, or labels from file; 2D or 3D.
- Foci with a seeded watershed on a filtered image (threshold in units of the image noise) or Spotiflow;
  colocalization between foci channels with a chance control; tracking of cells in time-lapses.
- Cells can be removed from the analysis by hand; tables, label images and images for Fiji with overlays.
- Reads TIFF, CZI and LIF. Installers for Windows (setup.exe), macOS and Linux, with the neural networks of
  your choice; updates from within the program.
- **Estimate** per foci channel: the foci size and threshold from the shown image, as a starting point.
- The viewer can show the foci channels after the filter (**Filtered**, key R), in percent of the threshold.
- *Cells to check* lists cells on the image edge, much larger or smaller than the others, or with far more foci;
  images with many saturated foci get a remark.
- The Z and T sliders have a handle as wide as one slice or frame. A browser tab whose program stopped (for
  example closed by the setup) says so, and comes back by itself when Foci Analyzer starts again.
- License: GNU General Public License v3.0 or later.
- Manual: https://bioimaging-nki.github.io/Foci-analyzer/

## Development before 2.0.0b1

First working version (0.1.0 in development). There are no earlier results to compare with.

- Viewer: WebGL, tiles per pyramid level, smooth zoom, slice or max projection, contrast and color per
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
  chance control in gray. The chance control is now on by default (10 repeats) and in both directions: the
  cell table gets `chance_...` columns. `Color cells by` fills the cells with a color for any of their numbers.
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
- Cells segmented in 3D: foci are detected in 3D too (projections are grayed out); before, a projection
  could put foci where no cell is in that slice.
- A foci size smaller than a pixel no longer stops the analysis (IndexError).
- Smooth outlines, cell numbers, color cells by label, an adjustable range for the color scale, the
  selected cell highlighted, the distance outside cells as a gray line.
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
  grayed out. The automatic Cellpose model works with Cellpose 3 and 4.
- Windows: `Foci Analyzer.exe` in the installation folder, with the logo as icon (also for the shortcuts);
  Linux: the icon is also put in the icon theme.
- **Files named** matches with or without the extension (`*.tif` found nothing before).
- A user manual (`manual/`, not online yet).
- Logo in the window, the browser tab, the setup and the shortcuts.
- Colocalization map: a white line with a black edge, and overlap hatched black and white, visible on any color.
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
  For Fiji, a green tick and Show in folder; *Color cells by* in groups; the Threshold help explains the noise.
- The progress bar fills up when a task ends and then empties at once (it slid back); it no longer jumps back
  during an analysis. The starting page no longer reloads itself all the time (the logo blinked).
- Settings page: installing a network no longer asks the other networks for their version while it runs (red
  errors) and hides its messages when it succeeded.
- Setup: networks taken from your own environment are remembered and shown as such.
- Fixed: before the first analysis (or with Live off) only the first change of a setting was kept; later
  changes went into an old copy and were lost. Removed cells can be selected with a Shift box again (to put
  them back). With Live on, another time frame of an image analyzed in the viewer is analyzed too.
- *Cells to check*: a cell is flagged for saturated foci only when more than a fifth of its foci are saturated.
  *Color cells by*: clearer names, and the channel in every foci entry.
