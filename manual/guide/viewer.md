# The viewer

The viewer shows the image with the result of the analysis on top. It loads large images in tiles, at the
resolution that the zoom needs.

![The Analyze page](../images/analyze.png)

## Moving around

- **Mouse wheel**: zoom; **drag**: move; **Fit** (++0++): the whole image.
- **Slice** and **Max projection** (++p++ switches): one Z slice or the brightest value over Z. When the foci
  are detected on another projection (extended depth of focus or sum), a third button shows that one.
- **Z** slider (++comma++ and ++period++): the slice. **T** slider (++bracket-left++ and ++bracket-right++):
  the time frame of a time-lapse; in an image without frames these keys step through the slices.
- The image list above the viewer chooses the image. A ✓ marks analyzed images: green, in the results
  folder with the current settings; blue, analyzed in the viewer with the current settings; grey, analyzed
  with other settings.

All keys: [Keys](../reference/keys.md).

## Channels

<figure markdown>
![The display panel](../images/display_panel.png){ width="280" }
</figure>

Per channel: the box shows or hides it (keys ++1++ to ++9++), the coloured square chooses its colour, and the
two ends of the slider set the darkest and brightest value shown. **Auto** sets the contrast from the image.
This changes the display only, never the measurements. The display stays the same when you open another
image with the same channels.

## Overlays

- **Cells** (++c++): the outlines of the nuclei or cells. When *Distance outside cells* is set, a thin
  dashed grey line shows how far outside a cell foci still count for it.
- **Foci** (++f++): the detected foci, lightly filled, in the colour of their foci channel; per channel
  they can be switched off.
- **Markers** (++m++): a dot on every focus. The dots follow the threshold slider at once; the outlines
  when the analysis has been repeated. When zoomed in, a ring with a black and white edge marks a focus that
  colocalizes with a focus of another channel.
- **Numbers** (++n++): the number of every cell, as in the tables, when it is large enough on the screen.
- **Smooth outlines** (++s++): outlines as smooth lines instead of along the pixel edges.
- **Colocalization** (++o++): colocalizing foci get a white outline with a black edge; where foci of two
  channels lie on top of each other is hatched.
- **Colour cells by**: fills every cell with a colour for one of its numbers (for example the number of
  foci), from dark blue (lowest) to yellow (highest). Drag the two ends of the slider below the colour
  scale to change its range. *label* gives every cell its own colour, which shows where one cell ends and
  the next begins; useful to judge a 3D segmentation slice by slice.

## The selected cell

Click a cell: its outline turns yellow, and **Selected cell** on the right lists its numbers (the same as
in `cells.tsv`). **This image** above it gives the totals of the image: cells, foci per channel, the noise
and threshold, and the colocalization.

<figure markdown>
![The numbers of the selected cell](../images/selected_cell.png){ width="280" }
</figure>

## Removing cells from the analysis

A cell that should not count (a mitotic or dead cell, two nuclei found as one, debris) can be left out:
select it and click **Remove from the analysis** (or press ++x++). It gets a red outline and leaves the
tables, together with its foci. Click it and press the button again (**Put back in the analysis**) to
undo.

Removed cells are remembered in the results folder, in `removed_cells.tsv`, as the points you clicked. So a
removed cell stays removed when the image is analyzed again, even with other segmentation settings: the
cell at that point is left out. Run **Analyze all images** again to update the tables in the results
folder. A different results folder starts without removed cells.

## Saving pictures

The **Save** section writes into the results folder:

- **View**: a PNG of the viewer as you see it (folder `pictures`).
- **Whole image**: the shown slice or projection at full resolution, with the channels and overlays as
  shown, as an RGB TIFF (folder `pictures`).
- **For Fiji**: the whole image (all channels, slices and frames) as a TIFF that opens in Fiji with the
  outlines of the cells and foci as overlay (folder `images`). It is saved in a `.zip`, about four times
  smaller; Fiji opens the `.zip` directly.
