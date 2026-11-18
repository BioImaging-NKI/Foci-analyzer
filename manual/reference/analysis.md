# How the analysis works

Every image (and every frame of a time-lapse) goes through the same steps. The viewer and **Analyze all
images** run the same code, so the preview shows exactly what the tables will hold (except with *Visible
part only*).

## 1. Preparation

*Remove border* pixels are cut from every edge, and with *XY binning* n × n pixels are averaged. The pixel
size and Z spacing come from the file; a stack without a Z spacing gets the pixel size, with a remark.

## 2. Nuclei or cells

The nuclei channel, or the cell channel when one is set, is segmented:

- **StarDist**, **Cellpose**, **InstanSeg**: neural networks, each in a Python environment of its own (see
  [The Settings page](../guide/settings-page.md)). Images sampled finer than 0.25 µm per pixel are made
  smaller for the network, and the outlines are scaled back. With a cell channel, Cellpose and InstanSeg
  also get the nuclei channel. StarDist works best on nuclei of a few tens of pixels: when its nuclei come
  out smaller than 16 or larger than 120 pixels, the image is also tried at ¼, ½, 2 and 4 times that size, and
  the result with nuclei closest to 40 pixels is kept (a scale far too large gives small fragments, so it
  is never chosen). With *Expected diameter* set, the image is scaled so nuclei are 40 pixels at once.
- **threshold**: the image is smoothed, Otsu's threshold (times the *Threshold* factor) separates
  objects from background, holes are filled, and touching objects are split by a watershed on the distance
  map.
- **labels from file**: label images made elsewhere, with the name of the image.

*Segment in* `2D` segments the maximum projection; the outlines hold for every slice. `3D` segments the
stack: Cellpose in real 3D for stacks of 20 slices or more; on thinner stacks, and always with StarDist and
InstanSeg, every slice is segmented and objects in neighboring slices are linked when they overlap by at
least half of the smaller one ("2.5D").

Objects smaller than *Smallest diameter* or larger than *Largest diameter* (diameter of a circle with the
same area) are removed, and with *Remove objects on the image edge* also objects that touch the edge.

## 3. Foci

Per foci channel, in the stack (*Detect foci in* `3D`) or in a projection (maximum, sum, or extended depth
of focus, which takes every pixel from the slice where it is sharpest):

1. **Filter.** The image is filtered at the *Foci size* (taken as the full width at half maximum of a
   focus; in Z times the *Axial elongation*): *DoG*, the difference of two Gaussians around that size;
   *LoG*, the Laplacian of Gaussian; or *wavelet*, the product of two neighboring à trous wavelet planes
   (Olivo-Marin, 2002). Spots of about the foci size stand out; smooth background disappears.
2. **Noise.** The noise is 1.4826 × the median absolute deviation of the filtered values inside the cells:
   a robust standard deviation that the foci themselves hardly affect. It is the typical size of the random
   ups and downs of the filtered image between the foci; for pure noise it equals the standard deviation. *Noise from* `per image` takes one
   value for all cells, `per cell` one per cell, and `fixed value` uses the same absolute threshold for all
   images.
3. **Candidates.** Local maxima of the filtered image inside the cells, with their border moved by *Cell
   border shift* (outwards when positive, inwards when negative), or anywhere with *Foci in the whole image*
   (then every focus counts for the nearest cell).
4. **Threshold.** A candidate is a focus when its height is at least *Threshold* × the noise. The viewer
   keeps all candidates, which is why the markers follow the threshold slider at once.
5. **Extent.** A seeded watershed grows every focus from its peak, down to half of its own peak height
   (*Focus extent* `half maximum`, so the size does not depend on the threshold) or to the threshold. A focus
   never grows into another focus or outside its cell. Foci outside *Smallest focus* and *Largest focus*
   are removed.

With **spotiflow**, the network Spotiflow finds the spots and gives each a probability; spots above
*Spotiflow probability* inside the cells are foci, and their extents are found as in step 5. A 2D Spotiflow
model on a stack finds the spots in the maximum projection and takes the slice of the brightest value.

Each focus is measured: position, size, mean and peak intensity, and the intensity summed above the local
background (the image without its foci: a gray-value opening with a window three focus diameters wide).

## 4. Colocalization

For every pair of foci channels, a focus of channel A colocalizes with channel B when it shares at least
*Smallest overlap* voxels with a B focus (`overlap`), or when the centers are closer than *Largest
distance* (`centroid distance`). Per cell this gives the number and share of colocalizing foci.

The **chance control** places the B foci at random positions inside the same cell, *Chance control
repeats* times, and gives the share that colocalizes then. Real colocalization is a share clearly above
that.

Intensity-based measures per cell: the Pearson correlation of the two channels, and Manders coefficients
(the share of the foci intensity of one channel inside foci of the other).

## 5. Time-lapse

With *Track cells* `overlap`, a cell is linked to the cell of the previous frame it overlaps most;
`trackastra` uses the network Trackastra, which also follows divisions. Every cell gets a `track` number.

## 6. Cells removed by hand

Cells removed in the viewer are left out of the tables, with their foci, after the analysis. They are
stored as points, so a new segmentation removes the cell that is at that point.
