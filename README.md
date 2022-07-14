# Foci-analyzer
ImageJ macro for the analysis of foci (e.g. DNA damage) in nuclei/cells

ImageJ macro to quantify foci in nuclei/cells. Works on 2D/3D images, including multiseries files, as long as all series have the same dimensions etc.
Foci can be measured in two channels, as well as colocalization.

More details on the workflow will follow soon.

![image](https://user-images.githubusercontent.com/68109112/179017837-4946a9ee-2602-4c52-a623-46dca8336c09.png)


## Installation / Requirements
► Requires the following Fiji update sites:
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


Author: Bram van den Broek, The Netherlands Cancer Institute, b.vd.broek@nki.nl
