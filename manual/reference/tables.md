# Table columns

Positions and sizes are in micrometres (µm) when the image is calibrated, otherwise in pixels.
Intensities are in the gray values of the image. *N*, *A* and *B* stand for channel numbers.

## cells.tsv: one row per nucleus or cell

| Column | Meaning |
|---|---|
| `condition`, `image`, `folder` | where the cell comes from |
| `frame` | the time frame (1 for a single image) |
| `cell` | the number of the cell in its image and frame, as shown by **Numbers** in the viewer |
| `x`, `y` | the center of the cell (µm from the top left corner of the analyzed image) |
| `area` | the area of the cell in the projection (µm²) |
| `volume` | the volume (µm³); only for cells segmented in 3D |
| `mean_chN` | mean intensity of channel N in the cell (in 2D segmentation: of the mean over the slices) |
| `foci_chN` | the number of foci of channel N |
| `foci_chN_total_intensity`, `foci_chN_mean_intensity` | the summed and mean integrated intensity of those foci, above the local background |
| `foci_chN_total_size`, `foci_chN_mean_size` | the summed and mean size of those foci (µm³ in 3D, µm² in 2D) |
| `pearson_chA_chB` | Pearson correlation of the two foci channels inside the cell (on the image the foci were detected in) |
| `coloc_chA_chB` | the number of A foci that colocalize with a B focus |
| `fraction_chA_with_chB` | the share of the A foci that colocalize with a B focus (empty when the cell has no A foci) |
| `chance_chA_with_chB` | the same share when the B foci are placed at random places in the cell: the chance control (*Chance control repeats*) |
| `manders_chA_in_chB` | Manders coefficient: the share of the foci intensity of A that lies inside B foci |
| `track`, `parent` | in a tracked time-lapse: the track of the cell, and the track it divided from |
| `check` | remarks that call for a look: *on edge* (the cell touches the image edge), *saturated foci in channel N* (more than a fifth of its foci reach the highest gray value; the counts are fine, their intensities are not) |
| `path` | the image file (and series) the row belongs to |

Cells removed by hand in the viewer are not in the table.

## foci.tsv: one row per focus

| Column | Meaning |
|---|---|
| `condition`, `image`, `folder`, `frame` | where the focus comes from |
| `cell` | the cell it belongs to (the `cell` of cells.tsv) |
| `channel`, `focus` | its foci channel, and its number in that channel |
| `track` | in a tracked time-lapse: the track of its cell |
| `x`, `y`, `z` | its center (µm; z 0 for foci found on a projection) |
| `size`, `voxels` | its size (µm³ in 3D, µm² in 2D) and number of pixels or voxels |
| `mean`, `peak` | its mean and highest intensity |
| `integrated` | its summed intensity above the local background |
| `height` | the height of its peak in the filtered image, in units of the threshold (× noise), or the Spotiflow probability |
| `with_chB` | whether it colocalizes with a focus of channel B |
