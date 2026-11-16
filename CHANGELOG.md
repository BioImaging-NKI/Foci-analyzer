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
- Status bar at the bottom: what the program is doing, step and progress bar, and when it is ready.
  Shortcut keys are shown next to their buttons; `[` and `]` step through slices when there is no time-lapse.
  Results do not change.
