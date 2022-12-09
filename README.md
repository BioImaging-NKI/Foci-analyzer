# Foci-analyzer

ImageJ macro for the analysis of foci (e.g. DNA damage) in nuclei (or cells). Works on 2D/3D fluorescence images, including multiseries files, as long as all series have the same dimensions.
The macro doesn't work on timelapse images (yet). A (slightly tedious) workaround is splitting your timepoints first. This option may be included in a future version.

Author: Bram van den Broek, The Netherlands Cancer Institute (b.vd.broek@nki.nl or bioimaging@nki.nl)

![image](https://user-images.githubusercontent.com/68109112/180581530-dd326026-cc74-4ce1-8d97-14518bfd4d73.png)

## Workflow summary
1. Nuclei are segmented (in a 2D projection) using the pre-trained deep learning network [StarDist](https://imagej.net/plugins/stardist). Alternatively, classic thresholding + watershedding can be used (though no parameters can be changed). As a third option, the deep learning network [Cellpose](https://github.com/MouseLand/cellpose) can be used to segment whole cells, thanks to the [Cellpose wrapper for Fiji](https://github.com/BIOP/ijl-utilities-wrappers) by BIOP.[^1]
[^1]: Currently, Cellpose is run using the 2D 'cyto' model on a single channel, with the default options and automatic diameter detection. Look for `run("Cellpose Advanced")` in the code and change parameters as seen fit.

2. Foci are detected in each nucleus in 2D or 3D. After Difference of Gaussians background subtraction, local maxima are detected and used as seeds for [MorpholibJ](https://imagej.net/plugins/morpholibj)'s [marker-controlled watershed](https://imagej.net/plugins/marker-controlled-watershed) (executed on the GPU using [CLIJ2/CLIJx](https://clij.github.io/)). Additionally, AreaMaxima local maxima detection can be used as detection method.
Thresholds are automatically calculated per image as 3 times the median of the standard deviations of the (outlier-removed) foci signal in the nuclei, and can be adapted using the threshold bias sliders.

3. Foci are quantified in a single channel A, or in two channels A and B.
For each nucleus a number of metrics are reported in a table: foci count, as well as mean/median values for foci intensity, area/volume, and whole nucleus intensity.
If two channels are selected, foci colocalization between channels is automatically calculated.
The table is saved as a `.tsv` file, which can be easily opened in other programs (e.g. Excel).

4. Segmented nuclei and foci are visualized as overlays on the original images for easy inspection. If desired, foci detection settings (threshold bias, min/max size) can be adapted before processing all input images. A colocalization map is also created when two channels are measured.

## Installation / Requirements
In order to run the macro, download the [latest release](https://github.com/BioImaging-NKI/Foci-analyzer/releases/download/v1.2/Foci_Analyzer_1_2.ijm) into a sensible folder, then drag&drop the `.ijm` file on the Fiji window. This will open the editor. Click 'Run' in the bottom left.

► [Add](https://imagej.net/update-sites/following) the following Fiji update sites:
- 3D ImageJ Suite
- CLIJ
- CLIJ2
- CLIJx-assistent
- CLIJx-assistent-extensions
- CSBDeep
- IJPB-plugins
- SCF MPI CBG
- StarDist

► If you want to also use Cellpose segmentation you also need:
- A working Cellpose 2.0 environment in Python
- PTBIOP update site, with proper settings. See https://github.com/BIOP/ijl-utilities-wrappers/blob/master/README.md#cellpose

## Short usage manual

The macro starts with a large dialog containing all options and parameters (click to enlarge):

<img src="https://user-images.githubusercontent.com/68109112/205070945-f84c1102-9121-437c-908e-3332dd3a620a.png" width="400">

The dialog has several sections, which are discussed below. All settings in this dialog will be remembered after you click `OK`.

### File handling and channel settings
- _Input files_ : Here you can specify which files to analyze, by adding them to the list, or drag&drop from a file explorer window
- _Output folder_ : the folder where all the analyzed images and results will be written.
- _Nuclei channel_ : The image channel on which nuclei segmentation will be performed. This segmentation always happens in 2D; for 3D images on a maximum intensity projection. (default: 1)
- _Foci channel A_ : the first foci channel (default: 2)
- _Foci channel B_ : the second foci channel (default: 3)
- _Also detect foci in Channel B and perform colocalization?_ : If checked, foci in both channels *A* and *B* will be analyzed, followed by a simple colocalization analysis. (default: checked)
- _Remove image borders (pixels)_ : Possibility to remove edges from the image. This can in particular be useful when the image edges have very different intensities, causing incorrect automatic nuclei segmentation.  (default: 0)
- _Image XY binning before analysis_ : Optional pixel binning in case the resolution is very high and the foci consist of many pixels. A value of 2 means: 2x2 pixels will be binned into 1 pixel. This reduces noise in the image and speeds up analysis. (default: 1)

### Nuclei segmentation settings
- _Nuclei segmentation method_ :
  - *Stardist nuclei segmentation* (default) uses the pretrained convolutional neural network [StarDist](https://github.com/stardist/stardist#readme) to recognize cell nuclei in fluoresence microscopy images. In general this method works very well on a large variety of samples.
  - *Cellpose cytoplasm segmentation* uses the deep learning network [Cellpose](https://github.com/MouseLand/cellpose#--cellpose--) (model: cyto) to recognize whole cells. Use this option if you want to measure foci in entire cells and/or you do not have a nuclear staining (but it can also work well for nuclei). Cellpose requires a few additional installations (see [Installation / Requirements](https://github.com/BioImaging-NKI/Foci-analyzer/edit/main/README.md#installation--requirements)).
  - _Classic nuclei segmentation_ allows the user to segment nuclei using manual/automatic thresholding is provided for legacy reasons. The method is almost always outperformed by the two other methods.
- _Stardist nuclei binning factor [1-n]_ : Stardist is trained on medium resolution images, and generally performs well on images with pixel sizes around 0.5 µm. For images with much smaller pixel size StarDist tends to 'oversegment' nuclei. In this case, increase the StarDist binning factor. (N.B. This option only affects the nuclei segmentation; it is different from the previously mentioned 'XY binning' parameter, which also changes the pixel size of the foci channels.) (default: 1)
- _Probability threshold [0.0-1.0] (StarDist/Cellpose)_ : Lower values will accept more nuclei; higher values will be more stringent.
- _Remove nuclei with diameter smaller than (units)_ : Objects smaller than circles having an area corresponding to this diameter will be removed. 'Units' depends on the image, and will almost always be 'µm', or otherwise 'pixels' in case the pxiel calibration values are missing. (default: 4)
- _Remove nuclei with diameter larger than (units)_ : Likewise, but for large objects. (default: 50)
- _Exclude nuclei on image edges_ : When checked, nuclei that touch the image edge will not be analyzed.
- _Manual nuclei removal_ : allows the user to erase ill-segmented nuclei before analysis. (default: No thanks) Options:
  - *No thanks* means no manual nuclei editing
  - *Manually remove nuclei* : Remove nuclei by leftclicking them in the image with the mouse. Editings will be saved to a small text file in the output folder.
  - *load previously saved removals (from output folder)* : If you have edited the segmented nuclei on this image before, it will load the previous nuclei removals from the file in the specified output folder. (Hence, if you change the output folder parameter this option will not work.)
  
<More info soon...>

