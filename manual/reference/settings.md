<!-- Made by tools/manual/gen_reference.py from config.py: do not edit by hand. -->

# All settings

Every analysis setting of the panel on the left of the **Analyze** page, with its default and the text of its **?**. Settings marked *advanced* are under **Advanced** in their section. Settings files (**Save settings…**, and `settings.toml` in every results folder) use the names in `code`, under the section name in brackets.

## Image

In a settings file: `[image]`.

| Setting | Default | Description |
|---|---|---|
| **Nuclei channel**<br>`nuclei_channel` | 1 | The channel with the nuclear stain (for example DAPI). Foci are counted per nucleus. 0 = no nuclei channel. At least 0. |
| **Cell channel**<br>`cell_channel` | 0 | A channel that shows whole cells (cytoplasm or membrane). When set, cells are segmented instead of nuclei and foci are counted per cell. 0 = none. With both channels at 0 the whole image is treated as one cell. At least 0. |
| **XY binning**<br>`xy_binning`<br>*advanced* | 1 | Average n × n pixels before the analysis. Use it for images sampled much finer than the foci size; it reduces noise and time. At least 1. |
| **Remove border**<br>`crop_border`<br>*advanced* | 0 px | Pixels removed from every image edge before the analysis, for example a dark stitching border. At least 0. |

## Nuclei / cells

In a settings file: `[segmentation]`.

| Setting | Default | Description |
|---|---|---|
| **Method**<br>`method` | stardist | How nuclei or cells are found. StarDist, Cellpose and InstanSeg are neural networks and need the models to be installed. StarDist: nuclei only, fast and robust. Cellpose: nuclei or cells, slower. InstanSeg: nuclei or cells, fast; made for tissue as well. All three also segment stacks (see Segment in). 'threshold' smooths the image, applies an Otsu threshold and splits touching objects on the distance map. 'labels from file' reads label images made elsewhere. Choices: `stardist`, `cellpose`, `instanseg`, `threshold`, `labels from file`. |
| **Segment in**<br>`dimension` | 2D | 2D segments the maximum projection of the stack; the outlines then hold for all slices. This is faster and often more robust, and can be combined with foci detection in 3D. 3D segments the stack itself: Cellpose in real 3D for stacks of 20 slices or more; StarDist, InstanSeg and Cellpose on thinner stacks segment every slice and link the slices where objects overlap ("2.5D"); threshold works in 3D directly. Choices: `2D`, `3D`. |
| **Smallest diameter**<br>`min_diameter` | 4 µm | Objects with a smaller diameter (of a circle with the same area) are removed. At least 0. |
| **Largest diameter**<br>`max_diameter` | 50 µm | Objects with a larger diameter are removed. 0 = no limit. At least 0. |
| **Remove objects on the image edge**<br>`exclude_on_edges`<br>*advanced* | on | A nucleus that is cut by the image edge has an incomplete foci count. |
| **Model**<br>`model` | (empty) | The pretrained model of the network, or a model of your own (the last choice). Cellpose 3 and Cellpose 4 models run in environments of their own; models of a network that is not installed are greyed out (install it on the Settings page). The list of models is on [Models](models.md). Shown when *Method* is `cellpose` or `instanseg` or `stardist`. |
| **Expected diameter**<br>`diameter`<br>*advanced* | 0 µm | The typical diameter of the nuclei or cells. 0 = automatic: Cellpose estimates it; StarDist is given the image at 0.25 µm per pixel and, when the nuclei then come out smaller than 16 or larger than 120 pixels, at the scale where they are closest to 40 pixels. Set it when the automatic choice is off. At least 0. Shown when *Method* is `cellpose` or `instanseg` or `stardist`. |
| **Probability threshold**<br>`probability`<br>*advanced* | 0.5 | Lower finds more and larger objects. Range: 0 to 1. Shown when *Method* is `cellpose` or `instanseg` or `stardist`. |
| **Threshold**<br>`threshold_factor`<br>*advanced* | 1 × Otsu | For 'threshold': the Otsu threshold of the image times this factor. Higher makes the objects smaller and splits more of them; lower makes them larger. Range: 0.2 to 3. Shown when *Method* is `threshold`. |
| **Labels folder**<br>`labels_folder` | (empty) | For 'labels from file': a folder with label images (TIFF) that have the same name as the images. Shown when *Method* is `labels from file`. |

## Foci

In a settings file: `[foci_general]`.

| Setting | Default | Description |
|---|---|---|
| **Detect foci in**<br>`dimension` | 3D | For stacks. 3D detects foci in the stack with filters that work in X, Y and Z. A projection detects them in one plane, which is often more reliable for thin samples and few slices. Ignored for 2D images. Choices: `3D`, `max projection`, `extended depth of focus`, `sum projection`. |
| **Distance outside cells**<br>`outside_distance`<br>*advanced* | 0 µm | Foci up to this distance outside a nucleus or cell still count for it. -1 = every focus counts for the nearest cell. At least -1. |
| **Axial elongation**<br>`z_elongation`<br>*advanced* | 2 × | How much longer a focus is in Z than in XY because of the microscope optics. Used for the filter size in Z (3D detection). Range: 0.5 to 10. |

## Foci channel (one per foci channel)

In a settings file: `[[foci]]`.

| Setting | Default | Description |
|---|---|---|
| **Channel**<br>`channel` | 2 | The channel with the foci. At least 1. |
| **Method**<br>`method` | watershed | watershed: the image is filtered so that spots of the foci size stand out, local maxima above the threshold become foci, and a seeded watershed gives every focus its extent. spotiflow: a neural network finds the spots and gives each a probability; the extents are found in the same way. Choices: `watershed`, `spotiflow`. |
| **Foci size**<br>`size` | 0.5 µm | Expected diameter of a focus (full width at half maximum). It sets the filter scale and the smallest distance between two foci. At least 0.05. |
| **Threshold**<br>`sensitivity` | 8 × noise | A local maximum is a focus when it is this many times higher than the noise of the filtered image. The noise is 1.4826 × the median absolute deviation of the filtered values inside the cells, which is hardly influenced by the foci themselves. Lower finds more foci. The viewer follows this slider directly. Range: 1 to 50. Shown when *Method* is `watershed`. |
| **Spotiflow probability**<br>`probability` | 0.5 | Spots below this probability are left out. Range: 0 to 1. Shown when *Method* is `spotiflow`. |
| **Filter**<br>`filter`<br>*advanced* | DoG | DoG: difference of two Gaussians around the foci size. LoG: Laplacian of Gaussian at the foci size; slightly sharper, more sensitive to noise. wavelet: product of two à trous wavelet planes; good for foci of mixed size on uneven background. All work in 3D for 3D detection. Choices: `DoG`, `LoG`, `wavelet`. Shown when *Method* is `watershed`. |
| **Focus extent**<br>`extent`<br>*advanced* | half maximum | Where a focus ends: at half of its own peak height (size does not depend on the threshold), or at the threshold. Choices: `half maximum`, `threshold`. Shown when *Method* is `watershed`. |
| **Smallest focus**<br>`min_size`<br>*advanced* | 2 voxels | Smaller foci are removed. At least 1. |
| **Largest focus**<br>`max_size`<br>*advanced* | 9999 voxels | Larger foci are removed. At least 1. |
| **Noise from**<br>`threshold_scope`<br>*advanced* | per image | per image: one noise value from all cells of the image. per cell: every cell its own noise. fixed value: the same absolute threshold for all images, which is the fair choice when conditions differ strongly in foci numbers. Choices: `per image`, `per cell`, `fixed value`. Shown when *Method* is `watershed`. |
| **Fixed threshold**<br>`fixed_threshold`<br>*advanced* | 0 | For 'fixed value': the threshold in units of the filtered image. The viewer shows the value that the current image uses. At least 0. Shown when *Method* is `watershed` and *Noise from* is `fixed value`. |
| **Spotiflow model**<br>`spotiflow_model` | general | The pretrained Spotiflow model, or a model of your own (the last choice: its folder). 3D models need a stack and 'Detect foci in' 3D. The list of models is on [Models](models.md). Shown when *Method* is `spotiflow`. |

## Colocalization

In a settings file: `[colocalization]`.

| Setting | Default | Description |
|---|---|---|
| **Colocalize foci channels**<br>`enabled` | on | Used when there are two or more foci channels; every pair of channels is compared. |
| **Method**<br>`method` | overlap | overlap: two foci colocalize when they share voxels. centroid distance: when their centres are closer than the largest distance. Choices: `overlap`, `centroid distance`. |
| **Smallest overlap**<br>`min_overlap`<br>*advanced* | 1 voxels | At least 1. Shown when *Method* is `overlap`. |
| **Largest distance**<br>`max_distance`<br>*advanced* | 0.3 µm | At least 0. Shown when *Method* is `centroid distance`. |
| **Chance control repeats**<br>`chance_repeats`<br>*advanced* | 10 | The colocalization expected by chance: the foci of the second channel are placed at random positions in the cell this many times (0 = off). The Results page shows it in grey next to the measured colocalization. At least 0. |

## Time-lapse

In a settings file: `[timelapse]`.

| Setting | Default | Description |
|---|---|---|
| **Track cells**<br>`tracking` | off | Follow nuclei or cells over the frames. overlap: a cell is linked to the cell it overlaps most in the previous frame. trackastra: a neural network that also handles divisions. Choices: `off`, `overlap`, `trackastra`. |
