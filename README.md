# Foci Analyzer (Python)

Detection, quantification and colocalization of foci in nuclei and cells, in a browser window on your own
computer. Python successor of the [Foci Analyzer](https://imagej.net/plugins/foci-analyzer) Fiji macro.
Under construction: `plan.md` has the design and the state per milestone.

## Install

No Python is needed beforehand; an internet connection is. Everything goes into one folder, and removing it
removes the program. Run the installer again to update or to add networks.

- **Windows**: run `FociAnalyzer-Setup-<version>.exe`. It asks which neural networks to install (StarDist,
  Cellpose, InstanSeg, Spotiflow, Trackastra; about 7 GB together) and whether you want to use networks from
  an environment you already have. Windows may warn about an unknown publisher: choose *More info* > *Run
  anyway*.
- **macOS and Linux**: unpack the program, open a Terminal in its folder and type `sh install.sh`. It asks the
  same questions, installs in `~/FociAnalyzer` and adds Foci Analyzer to Applications (macOS) or to the
  applications menu (Linux).

To check an installation: `<install folder>/.venv/Scripts/python.exe installer/check.py` (Windows; on macOS
and Linux `.venv/bin/python`) runs every installed network on a made-up image.

## Use

1. **Experiment**: open images or a folder, or make one condition per group you want to compare.
2. **Analyze**: choose the channels, press `Analyze this image` and move the threshold slider while you look
   at the image. The `?` next to every setting explains it. The panel on the right has the channels, the
   overlays (cells, foci, markers, the colocalization map, a colour per cell for any of its numbers) and the
   numbers of the image and of a clicked cell. `Analyze all images` writes the tables to the results folder.
3. **Results**: graphs per condition with tests (also colocalization against chance, and tracks of a
   time-lapse), a list per image and the cells to check. A click on a point opens that cell in the viewer.
4. **Settings**: which networks are installed, whether they use the graphics card, and which Python runs
   them (your own Cellpose environment can be given there).

Images: TIFF (ImageJ, OME), Zeiss CZI and Leica LIF; 2D, stacks, time-lapse.

## For developers

```
uv sync
uv run foci-analyzer
uv run pytest -q
```

Installers: `installer/install.py` does the work on every system; `installer/build-windows.ps1` builds the
Windows setup (needs Inno Setup 6). The GitHub workflow `build` tests on Windows, Linux and macOS, and
builds and tries the installers when run by hand or for a version tag.

The networks live in separate environments next to the program environment: `models-env` (PyTorch: Cellpose,
InstanSeg, Spotiflow, Trackastra) and `stardist-env` (TensorFlow 2.15: StarDist). Another Python can be given
on the Settings page, or with `FOCI_ANALYZER_MODELS_PYTHON`.
