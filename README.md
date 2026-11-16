# Foci Analyzer (Python)

Detection, quantification and colocalization of foci in nuclei and cells, in a browser window on your own
computer. Python successor of the [Foci Analyzer](https://imagej.net/plugins/foci-analyzer) Fiji macro.
Under construction: `plan.md` has the design and the state per milestone.

## Install (Windows)

Double-click `install.bat`. It needs an internet connection and no Python: it fetches
[uv](https://docs.astral.sh/uv/), installs the program in `%LOCALAPPDATA%\FociAnalyzer`, asks whether to add
the neural networks (Cellpose, InstanSeg, Spotiflow, Trackastra; about 6 GB) and puts a shortcut on the
desktop. Run it again to update or to add the networks later.

## Use

1. **Experiment**: open images or a folder, or make one condition per group you want to compare.
2. **Analyze**: choose the channels and move the threshold slider while you look at the image. The `?` next
   to every setting explains it. `Analyze all images` writes the tables to the results folder.
3. **Results**: graphs per condition with tests, a list per image and the cells to check. A click on a point
   opens that cell in the viewer.

Images: TIFF (ImageJ, OME), Zeiss CZI and Leica LIF; 2D, stacks, time-lapse.

## For developers

```
uv sync
uv run foci-analyzer
uv run pytest -q
```

The networks live in a separate environment (`models-env` next to the program environment, or the path in
`FOCI_ANALYZER_MODELS_PYTHON`, or the setting *Python of the models*).
