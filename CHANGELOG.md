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
- Status bar at the bottom: what the program is doing, step and progress bar, and when it is ready.
  Shortcut keys are shown next to their buttons; `[` and `]` step through slices when there is no time-lapse.
  Results do not change.
