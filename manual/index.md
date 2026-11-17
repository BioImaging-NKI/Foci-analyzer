# Foci Analyzer

Foci Analyzer detects and quantifies foci (for example DNA-damage foci of γH2AX, 53BP1 or RPA) in nuclei or
whole cells of 2D and 3D fluorescence microscopy images. It segments the nuclei or cells, detects the foci
of one or more channels, measures their number, size and intensity per cell, and determines the
colocalization between foci channels. Conditions are compared on the Results page, with tests that take
into account that cells of one image are not independent.

It is the Python successor of the Fiji macro
[Foci Analyzer](https://github.com/BioImaging-NKI/Foci-analyzer) and runs in your web browser, on your own
computer: images and results never leave it.

![The Analyze page: settings on the left, the image with the detected nuclei and foci in the middle, display and results on the right](images/analyze.png)

## What it does

- **Nuclei or cells**: StarDist, Cellpose (3 and 4), InstanSeg, a threshold, or label images made
  elsewhere; in 2D, or in 3D for stacks.
- **Foci**: a seeded watershed on a filtered image (difference of Gaussians, Laplacian of Gaussian or
  wavelets), with a threshold in units of the image noise, or the neural network Spotiflow; in 3D or on a
  projection.
- **Colocalization** of the foci of every pair of channels, by overlap or distance, with a chance control,
  and intensity-based Pearson and Manders coefficients per cell.
- **A live preview**: the result appears in the image while you change the settings; the foci markers
  follow the threshold slider at once.
- **Experiments**: conditions with replicate folders, all images analyzed side by side, tables per cell and
  per focus, graphs and statistics.
- **Time-lapses**: cells tracked over time by overlap or with Trackastra.
- Reads TIFF, CZI (with scenes) and LIF (with series); images of any size, also large stitched ones.

## Where to start

- [Installation](start/install.md): Windows, macOS and Linux.
- [Your first analysis](start/first-analysis.md): a walk-through from images to a graph.
- [All settings](reference/settings.md): every setting with its default and explanation.
- [How the analysis works](reference/analysis.md): the methods.
