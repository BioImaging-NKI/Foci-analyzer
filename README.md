# Foci-analyzer

ImageJ macro for the analysis of foci (e.g. DNA damage) in nuclei (or cells). Works on 2D/3D fluorescence images, including multiseries files, as long as all series have the same dimensions.

## Quick workflow
1. Nuclei are segmented using the pre-trained deep learning network [StarDist](https://imagej.net/plugins/stardist). Alternatively, classic thresholding + watershedding can be used (though no parameters can be changed). As a third option, the deep learning network [Cellpose](https://github.com/MouseLand/cellpose) can be used to segment whole cells, thanks to the [Cellpose wrapper for Fiji](https://github.com/BIOP/ijl-utilities-wrappers) by BIOP.[^1]
[^1]: Currently, Cellpose is run using the 'cyto' model on a single channel, with the default options. Look for `run("Cellpose Advanced")` in the code and change parameters as seen fit.

2. Foci are quantified in a single channel A, or in two channels A and B.
For each nucleus the foci count is reported, as well as mean/median values for foci intensity, area/volume, and whole nucleus intensity.
If two channels are selected With the latter, foci colocalization between channels is automatically calculated.

N.B. The macro doesn't work on timelapse images (yet). A (slightly tedious) workaround is splitting your timepoints first. This option may be included in a future version.

![image](https://user-images.githubusercontent.com/68109112/179017837-4946a9ee-2602-4c52-a623-46dca8336c09.png)

## Installation / Requirements
In order to run the macro, download the [latest release](https://github.com/BioImaging-NKI/Foci-analyzer/releases/tag/v1.0), then drag&drop `Foci_analyzer.ijm` on the Fiji window. This will open the editor. Click 'Run' in the bottom left.


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
- A working Cellpose Python environment
- PTBIOP update site, with proper settings. See https://github.com/BIOP/ijl-utilities-wrappers/blob/master/README.md#cellpose

## Short manual for using the macro in Fiji

To run the macro, drag&drop it on the Fiji window. This will open the editor. Click 'Run' in the bottom left.
A large dialog with options and parameters will appear (click to enlarge):

<img src="https://user-images.githubusercontent.com/68109112/180569141-d6b79331-8ee5-4561-b9a3-0e1a7e9ba659.png" width="400">


Author: Bram van den Broek, The Netherlands Cancer Institute, b.vd.broek@nki.nl


# Markdown Cheat Sheet

Thanks for visiting [The Markdown Guide](https://www.markdownguide.org)!

This Markdown cheat sheet provides a quick overview of all the Markdown syntax elements. It can’t cover every edge case, so if you need more information about any of these elements, refer to the reference guides for [basic syntax](https://www.markdownguide.org/basic-syntax) and [extended syntax](https://www.markdownguide.org/extended-syntax).

## Basic Syntax

These are the elements outlined in John Gruber’s original design document. All Markdown applications support these elements.

### Heading

# H1
## H2
### H3

### Bold

**bold text**

### Italic

*italicized text*

### Blockquote

> blockquote

### Ordered List

1. First item
2. Second item
3. Third item

### Unordered List

- First item
- Second item
- Third item

### Code

`code`

### Horizontal Rule

---

### Link

[Markdown Guide](https://www.markdownguide.org)

### Image

![alt text](https://www.markdownguide.org/assets/images/tux.png)

## Extended Syntax

These elements extend the basic syntax by adding additional features. Not all Markdown applications support these elements.

### Table

| Syntax | Description |
| ----------- | ----------- |
| Header | Title |
| Paragraph | Text |

### Fenced Code Block

```
{
  "firstName": "John",
  "lastName": "Smith",
  "age": 25
}
```

### Footnote

Here's a sentence with a footnote. [^1]

[^1]: This is the footnote.

### Heading ID

### My Great Heading {#custom-id}

### Definition List

term
: definition

### Strikethrough

~~The world is flat.~~

### Task List

- [x] Write the press release
- [ ] Update the website
- [ ] Contact the media

### Emoji

That is so funny! :joy:

(See also [Copying and Pasting Emoji](https://www.markdownguide.org/extended-syntax/#copying-and-pasting-emoji))

### Highlight

I need to highlight these ==very important words==.

### Subscript

H~2~O

### Superscript

X^2^
