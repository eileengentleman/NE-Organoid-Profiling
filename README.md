This repository contains code to analyse fluorescence intensity in elongated neuroepithelial organoids.


# Nuclear Beta-Catenin Intensity Quantification

Segments the nuclear region from the mid-plane of a C2 (nuclei) TIFF z-stack, then measures C1 (beta-catenin) intensity inside that mask across all z-planes.

## System Requirements

- Python 3.10+
- `numpy`
- `pandas`
- `pillow`
- `scipy`
- `scikit-image`
- No special hardware required; runs on any desktop OS.

## Installation

```bash
pip install numpy pandas pillow scipy scikit-image
````

Typical installation time: ~2–5 minutes on a normal desktop.

## Demo

### Input

An `<experiment-dir>/raw/` directory containing paired TIFF files named:

```text
C1-<name>(series N).tif
C2-<name>(series N).tif
```

The paired C1 and C2 TIFF stacks must have the same shape.

### Run

```bash
python Supplementary_Code_2.py \
    --experiment-dir /path/to/experiment \
    --output-dir /path/to/results
```

### Output

Results are written to `<output-dir>/`.

* `per_tiff_layer_csv/*_stack_intensity.csv`
  Contains measurements for each z-plane:

  * `stack_number`
  * `intensity`
  * `nuclei_number`
  * `target_value` — mean intensity per pixel

* `qc_overlays/*_qc.png`
  Nuclei image with the segmented mask boundary overlaid for quality control.

### Run Time

Seconds per series; typically a few minutes for a full dataset.

## Instructions for Use

### Arguments

| Argument               | Default        | Description                                             |
| ---------------------- | -------------- | ------------------------------------------------------- |
| `--experiment-dir`     | Script default | Folder containing the `raw/` subfolder                  |
| `--output-dir`         | Script default | Directory where results are written                     |
| `--min-area`           | `80`           | Minimum nucleus size in pixels                          |
| `--gaussian-sigma`     | `1.4`          | Gaussian smoothing applied before thresholding          |
| `--local-block-size`   | `151`          | Block size used for local thresholding                  |
| `--adaptive-threshold` | Off            | Enables combined local and global adaptive thresholding |
| `--start-series`       | None           | First series to process                                 |
| `--end-series`         | None           | Last series to process                                  |

## Segmentation Method

The nuclear mask is generated from the middle z-plane of the C2 stack using the following steps:

1. Apply Gaussian smoothing.
2. Apply Otsu thresholding.
3. Optionally combine Otsu thresholding with local adaptive thresholding when `--adaptive-threshold` is enabled.
4. Remove small objects.
5. Fill holes and remove small holes.
6. Apply the resulting nuclear mask from the middle C2 z-plane to every z-plane of the paired C1 stack.
7. Quantify beta-catenin intensity within the nuclear mask for each C1 z-plane.

