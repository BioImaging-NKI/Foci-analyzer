# Foci-analyzer

ImageJ macro for the analysis of foci (e.g. DNA damage) in nuclei (or cells). Works on 2D/3D fluorescence images, including multiseries files, as long as all series have the same dimensions.
The macro doesn't work on timelapse images (yet). A (slightly tedious) workaround is splitting your timepoints first. This option may be included in a future version.

Author: Bram van den Broek, The Netherlands Cancer Institute, b.vd.broek@nki.nl or bioimaging@nki.nl

![image](https://user-images.githubusercontent.com/68109112/180581530-dd326026-cc74-4ce1-8d97-14518bfd4d73.png)

## Workflow summary
1. Nuclei are segmented (in a 2D projection) using the pre-trained deep learning network [StarDist](https://imagej.net/plugins/stardist). Alternatively, classic thresholding + watershedding can be used (though no parameters can be changed). As a third option, the deep learning network [Cellpose](https://github.com/MouseLand/cellpose) can be used to segment whole cells, thanks to the [Cellpose wrapper for Fiji](https://github.com/BIOP/ijl-utilities-wrappers) by BIOP.[^1]
[^1]: Currently, Cellpose is run using the 2D 'cyto' model on a single channel, with the default options and automatic diameter detection. Look for `run("Cellpose Advanced")` in the code and change parameters as seen fit.

2. Foci are detected in each nucleus in 2D or 3D. After Difference of Gaussians background subtraction, local maxima are detected and used as seeds for [MorpholibJ](https://imagej.net/plugins/morpholibj)'s [marker-controlled watershed](https://imagej.net/plugins/marker-controlled-watershed) (executed on the GPU using [CLIJ2/CLIJx](https://clij.github.io/)). Additionally, AreaMaxima local maxima detection can be used as detection method.
Thresholds are automatically calculated as 3 times the median of the standard deviations of the (outlier-removed) nuclei, and can be adapted using the threshold bias parameters.

3. Foci are quantified in a single channel A, or in two channels A and B.
For each nucleus a number of metrics are reported in a table: foci count, as well as mean/median values for foci intensity, area/volume, and whole nucleus intensity.
If two channels are selected, foci colocalization between channels is automatically calculated.
The table is saved as a `.tsv` file, which can be easily opened in other programs (e.g. Excel).

4. Segmented nuclei and foci are visualized as overlays on the original images for easy inspection. If desired, foci detection settings (threshold bias, min/max size) can be adapted before processing all input images. A colocalization map is also created when two channels are measured.

## Installation / Requirements
In order to run the macro, download the [latest release](https://github.com/BioImaging-NKI/Foci-analyzer/releases/download/v1.1/Foci_Analyzer_1_1.ijm) into a sensible folder, then drag&drop the `.ijm` file on the Fiji window. This will open the editor. Click 'Run' in the bottom left.

► Requires the following [Fiji update sites](https://imagej.net/update-sites/following) to be installed:
- 3D ImageJ Suite
- CLIJ
- CLIJ2
- CLIJx-assistent
- CLIJx-assistent-extensions
- CSBDeep
- IJPB-plugins
- SCF MPI CBG
- StarDist

► In order to run Cellpose segmentation you also need:
- A working Cellpose 2.0 environment in Python
- PTBIOP update site, with proper settings. See https://github.com/BIOP/ijl-utilities-wrappers/blob/master/README.md#cellpose

## Short manual for using the macro in Fiji

The macro starts with a large dialog with options and parameters (click to enlarge):

<img src="https://user-images.githubusercontent.com/68109112/180569141-d6b79331-8ee5-4561-b9a3-0e1a7e9ba659.png" width="400">

...more info will follow soon...
