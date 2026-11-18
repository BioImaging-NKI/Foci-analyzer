# The results folder

**Analyze all images** writes into the results folder (chosen on the [Experiment page](experiment.md);
by default `Foci_Analyzer_results` next to the first image). A new run replaces the tables of the previous
one. Files are written under a temporary name and renamed when complete, so an interrupted run never
leaves half a table.

| File or folder | What it holds |
|---|---|
| `cells.tsv` | one row per nucleus or cell of every image and frame; see [Table columns](../reference/tables.md) |
| `foci.tsv` | one row per focus |
| `settings.toml` | the settings of this run; **Load settings…** reads it back |
| `run.json` | when the run was made, with which version, the images and any errors |
| `removed_cells.tsv` | the cells removed by hand in the viewer, as clicked points (image, frame, x, y, z in pixels; z = −1 for a click in a projection) |
| `labels/` | with **Also save label images**: per image `<condition>_<image>_cells.tif` and `<condition>_<image>_foci_ch<N>.tif` (T, Z, Y, X, 16-bit) |
| `images/` | with **Also save images with overlays (Fiji)**, or **For Fiji**: the image as an ImageJ TIFF with the outlines as overlay, in a `.zip` that Fiji opens directly |
| `pictures/` | **View** (PNG) and **Full image** (RGB TIFF) from the Save images section of the viewer |
| `graph_<measure>.tsv` | **Save graph data** on the Results page |
| `.viewer/` | the outlines and foci of the last run, so the viewer can show an image without analyzing it again; can be deleted |

The `.tsv` tables are tab-separated text. They open in Excel (double-click, or *Data → From Text*), in R
with `read.delim("cells.tsv")` and in Python with `pandas.read_csv("cells.tsv", sep="\t")`.

## Label images

The label images are in the grid of the analysis: after *XY binning* and without the *Remove border*
pixels (the pixel size in the file says so). Every cell has its own number, the same as in `cells.tsv`. In
a tracked time-lapse a cell keeps one number over time: the number of its track. Foci are numbered per
frame, as in `foci.tsv`. They open in Fiji and napari.

## Images for Fiji

Each image becomes a hyperstack with all its channels, slices and frames, the channel colours and contrast,
and the outlines as an overlay: cells in white, foci in the colour of their channel, each on its own slice
and frame (cells found in 2D are shown in every slice). Each is a `.zip` holding one TIFF: open the `.zip`
in Fiji directly (drag it onto Fiji, or *File → Open*). Zipped, the file is about four times smaller; a
compressed TIFF would be small too, but Fiji drops the overlay of a compressed TIFF.
