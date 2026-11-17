# Installation

Foci Analyzer installs into one folder of its own, with its own Python and, if you choose them, the neural
networks. Removing that folder removes everything. No administrator rights are needed.

The neural networks are optional: without them Foci Analyzer segments nuclei with a threshold and detects
foci with the watershed method. With them (about 7 GB, 10 GB with Cellpose 4) it finds nuclei and cells far
more reliably. They use an NVIDIA graphics card when there is one, the graphics chip of an Apple Silicon
Mac, or else the processor.

## Windows

1. Download `FociAnalyzer-Setup.exe` and run it.
2. Windows may warn that the publisher is unknown (the setup is not signed): choose **More info** and
   **Run anyway**.
3. Choose the networks. StarDist is recommended for nuclei; Cellpose 4 is not selected by default
   (large, and slow without a graphics card).
4. If you already have a Python environment with one of the networks (for example your own Cellpose), give
   its `python.exe` on the next page: that network is then not installed again.
5. The setup downloads Python and the networks; this takes 5 to 20 minutes.

Start Foci Analyzer from the start menu, or with `Foci Analyzer.exe` in the installation folder
(`%LOCALAPPDATA%\FociAnalyzer`). It opens in your web browser.

Silent installation, for example for several computers:
`FociAnalyzer-Setup.exe /SILENT /NETWORKS=stardist,spotiflow`
(`all`, `none`, or a list of `stardist`, `cellpose`, `cellpose4`, `instanseg`, `spotiflow`, `trackastra`).

## macOS and Linux

In a terminal, in the folder with the downloaded program:

```bash
sh install.sh
```

The script fetches [uv](https://docs.astral.sh/uv/) (a small tool that installs Python and packages),
asks which networks to install, and adds Foci Analyzer to Applications (macOS) or the applications menu
(Linux). No Apple developer account is involved: files fetched in the terminal are not blocked by macOS.
On a Mac with an Intel processor, CZI files cannot be read.

## Adding networks later

On the [Settings page](../guide/settings-page.md), a network that is not installed has an **Install**
button. You can also point Foci Analyzer to a Python environment of your own there.

## Removing

Windows: *Settings → Apps → Foci Analyzer → Uninstall*. macOS and Linux: remove the installation folder
(`~/FociAnalyzer` unless you chose another), the app in `~/Applications` or the menu entry, and
`~/.foci_analyzer` (the remembered settings).
