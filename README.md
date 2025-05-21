# C-Factor BRDF Correction in GEE w/ Python

## Overview

This repository provides a Python implementation of the c-factor BRDF logic in Google Earth Engine. The algorithm corrects Landsat reflectance data to nadir BRDF-adjusted reflectance (NBAR) using fixed MODIS-derived BRDF spectral model parameters. It accounts for viewing geometry, solar geometry, and band-specific BRDF coefficients to improve reflectance consistency across different Landsat scenes.

The Python implementation is based on a JavaScript implementation that was created at Daniel Wiell and Erik Lindquist of the United Nations Food and Agriculture Organization. Their implementation is available at this [link](https://code.earthengine.google.com/0bf07da7cdab0d0ae90962e9259ce8ec).

The algorithm is described in detail in Roy et al. 2016 ([paper link](https://www.sciencedirect.com/science/article/pii/S0034425716300220)): 

Roy, D. P., Zhang, H. K., Ju, J., Gomez-Dans, J. L., Lewis, P. E., Schaaf, C. B., ... & Kovalskyy, V. (2016). A general method to normalize Landsat reflectance data to nadir BRDF adjusted reflectance. Remote Sensing of Environment, 176, 255-271.

One modification of note is in the `compute_brdf` function. The original JavaScript implementation muliplied the term "kvol" by 3. If this not done than the effects of the BRDF correction are not apparent. In practice, this heuristic seems to work well and eliminates cross-track illumination effects.

## Usage

The core functionality for the package is the `apply_cfactor_brdf_correction` function. This function is designed to be mapped over a collection of Landsat satellite imagery. This function expects that the Landsat scenes have the TM/ETM+ naming scheme (i.e., "B1", "B2", "B3", "B4", "B5", and "B7"). If this is not desireable, you can easily modify the keys in the *coeffs_by_band* dictionary contained in the `apply_cfactor_brdf_correction` function.

## Installation

### A. The simple option

Download the repository and copy the `cfactor.py` script wherever you need it in your local code base or project. 

### B. Install as a module

1. Clone or download this repository.
2. Navigate to the repository.
3. Install it as a Python package:
   ```bash
   pip install .
4. Then, import the module in your scripts and use the `apply_cfactor_brdf_correction` function as needed.

## Implementation validation

To verify that this Python implementation produces the same results as as the JavaScript implementation a small test was conducted 

First, I selected a Landsat 8 scene in Brazil. This image was fed through the original JavaScript implementation of the c-factor BRDF algorithm. The GEE script used is available at this [link](https://code.earthengine.google.com/dadb5bffce05dc52f282b9e0688acd79).

Second, the same scene was fed through the Python implementation using the demonstration script `examples/usage_example.py`. 

The original image and the two BRDF-adjusted images can be visualized in the GEE playground using this [link](https://code.earthengine.google.com/07f30bed958fea3ed43446993a699200). If you use the inspector tool, you can see that at each pixel, the values of each band are the same for both images. Similarly, the pixel-level percent difference is 0 everywhere... which is a good sign!

## Improvements? Bugs?

If you have any suggestions for changes or if you identify a bug, please create a GitHub issue.
