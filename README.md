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
