# Installation

Foci Analyzer installs into one folder of its own, with its own Python and, if you choose them, the neural
networks. Removing that folder removes everything. No administrator rights are needed.

The neural networks are optional: without them Foci Analyzer segments nuclei with a threshold and detects
foci with the watershed method. With them it finds nuclei and cells far more reliably. Cellpose, InstanSeg,
Spotiflow and Trackastra share one PyTorch (4.4 GB with an NVIDIA graphics card, 1 GB without); StarDist
takes 1.6 GB and Cellpose 4, with a PyTorch of its own, 4.5 GB. The setup shows the total for your choice.
The PyTorch networks use an NVIDIA graphics card when there is one, or the graphics chip of an Apple Silicon
Mac; StarDist runs on the processor.

## Windows

1. Download `FociAnalyzer-Setup.exe` and run it.
2. Windows may warn that the publisher is unknown (the setup is not signed): choose **More info** and
   **Run anyway**.
3. Choose the networks. StarDist is recommended for nuclei; Cellpose 4 is not selected by default
   (large, and slow without a modern graphics card). The size of each and the total are shown.
4. If you already have a Python environment with one of the networks you left unticked (for example your
   own Cellpose), give its `python.exe` on the next page: Foci Analyzer then uses it. A venv and a conda
   environment (`...\envs\name\python.exe`) both work. This page only appears when a network is left
   unticked.
5. The setup downloads Python and the networks, and shows what it is doing; this takes 5 to 20 minutes.
   Everything is also written to `install-output.log` in the installation folder. A Foci Analyzer that is
   still running from that folder is closed first (your experiment and settings are kept).

Start Foci Analyzer from the start menu, or with `Foci Analyzer.exe` in the installation folder
(`%LOCALAPPDATA%\FociAnalyzer`). It opens in your web browser: first a page with the logo while the program
starts, then the window. Starting it a second time shows the window that is already open.

Running the setup when Foci Analyzer is installed already offers **Update it, or add networks** (into the
same folder; installed networks are ticked and stay) or **Remove it**. Networks that Foci Analyzer takes
from an environment of your own are marked *now from your own environment*; they stay unticked, and the next
page shows their `python.exe` again. Ticking one installs it and uses that instead of your own. Updating can also be done from Foci
Analyzer itself: see [Updates](../guide/settings-page.md#updates).

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
