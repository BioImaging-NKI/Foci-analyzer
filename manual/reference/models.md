<!-- Made by tools/manual/gen_reference.py from models.py: do not edit by hand. -->

# Models

The pretrained models that the **Model** lists offer, per network. A model needs its network to be installed (see [The Settings page](../guide/settings-page.md)); models of a network that is not installed are greyed out. *Own model…*, at the bottom of every list, takes the folder or file of a model you trained yourself.

## StarDist

More about these models: <https://github.com/stardist/stardist#pretrained-models-for-2d>

| Name | In the list | Runs in |
|---|---|---|
| `(default)` | Versatile fluorescent nuclei (2D_versatile_fluo) | StarDist |

## Cellpose

More about these models: <https://cellpose.readthedocs.io/en/latest/models.html>

| Name | In the list | Runs in |
|---|---|---|
| `(default)` | Automatic: nuclei or cyto3 (Cellpose 3), CellposeSAM (Cellpose 4) | Cellpose 3 |
| `nuclei` | Nuclei (Cellpose 3) | cellpose3 |
| `cyto3` | Cells, general (Cellpose 3) | cellpose3 |
| `tissuenet_cp3` | Tissue, cells and nuclei (Cellpose 3) | cellpose3 |
| `cpsam_v2` | CellposeSAM v2: any cell or nucleus (Cellpose 4) | Cellpose 4 |
| `cpdino-vitb` | CellposeDINO ViT-B (Cellpose 4; needs the dinov3 package, see the help) | Cellpose 4 |

## InstanSeg

More about these models: <https://github.com/instanseg/instanseg>

| Name | In the list | Runs in |
|---|---|---|
| `(default)` | Fluorescent nuclei and cells (fluorescence_nuclei_and_cells) | InstanSeg |
| `brightfield_nuclei` | Nuclei in brightfield (H&E, DAB) | InstanSeg |

## Spotiflow

More about these models: <https://weigertlab.org/spotiflow/pretrained.html>

| Name | In the list | Runs in |
|---|---|---|
| `general` | General: spots from many microscopes (2D) | Spotiflow |
| `hybiss` | HybISS spots (2D) | Spotiflow |
| `synth_complex` | Synthetic spots with aberrations and background (2D) | Spotiflow |
| `fluo_live` | Live-cell fluorescence (2D) | Spotiflow |
| `synth_3d` | Synthetic 3D spots (3D stacks) | Spotiflow |
| `smfish_3d` | smFISH in 3D (3D stacks) | Spotiflow |
