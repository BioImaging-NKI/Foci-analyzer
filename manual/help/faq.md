# Questions and problems

## Foci Analyzer does not open in the browser

Open `http://127.0.0.1:8770/` yourself. Foci Analyzer keeps running when you close the browser tab; **Quit**
at the top right stops it. Started a second time, it shows the window that is already running.

## The first analysis takes long

Starting a neural network takes up to a minute the first time (the status bar says so); the next analyses
take seconds. Foci Analyzer starts the networks in use when the window opens, to save that time.

## Are the networks using the graphics card?

The [Settings page](../guide/settings-page.md) says per network whether it runs on the graphics card or the
processor. StarDist always runs on the processor on Windows (fast enough; the version of TensorFlow that can
use the graphics card is a hundred times slower there).

## The nuclei are not found well

- Check the *Nuclei channel*.
- Try another method or model: StarDist is robust for nuclei; Cellpose and InstanSeg for cells and tissue.
- *Smallest diameter* and *Largest diameter* remove objects outside that range: check that your nuclei fit.
- For Cellpose, give the *Expected diameter* (under **Advanced**) when the automatic one is off.
- Use **Colour cells by** *label* to see where one object ends and the next begins.
- Single wrong nuclei can be [removed from the analysis](../guide/viewer.md#removing-cells-from-the-analysis).

## Too many or too few foci

- The *Foci size* matters as much as the *Threshold*: it sets the filter scale and the smallest distance
  between two foci. Measure a few foci (their width at half maximum) and use that.
- Move the *Threshold* slider and watch the markers: they follow at once.
- For conditions that differ strongly in foci numbers, *Noise from* `fixed value` gives the same absolute
  threshold for all images.

## Why do some markers have a ring?

A ring (with a black and white edge, visible when zoomed in) marks a focus that colocalizes with a focus of
another foci channel.

## The Fiji images have no overlay

The overlay is lost when the TIFF is saved compressed. Use the files as Foci Analyzer writes them.

## The results do not match the Fiji macro

Foci Analyzer uses its own methods (a threshold in units of the noise, foci extents at half maximum, 3D
detection); its counts are not meant to equal those of the Fiji macro.

## Where are my settings and removed cells?

The analysis settings are in `settings.toml` in every results folder, and remembered by the window.
Removed cells are in `removed_cells.tsv` in the results folder: with another results folder, no cells are
removed.
