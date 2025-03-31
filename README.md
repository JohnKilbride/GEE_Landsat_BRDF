# C-Factor BRDF Correction

## Overview

This repository provides a Python implementation of the c-factor BRDF logic. The approach corrects Landsat reflectance data to nadir BRDF-adjusted reflectance (NBAR) using fixed MODIS-derived BRDF spectral model parameters. 

It accounts for viewing geometry, solar geometry, and band-specific BRDF coefficients to improve reflectance consistency across different Landsat scenes.

The original JavaScript implementation was created at Daniel Wiell and Erik Lindquist of the United Nations Food and Agriculture Organization. The original implementation is available at this [link](https://code.earthengine.google.com/0bf07da7cdab0d0ae90962e9259ce8ec).

## Implementation validation

To verify that this Python implementation produces the same results as as the JavaScript implementation a small test was conducted 

First, I selected a Landsat 8 scene in Brazil. This image was fed through the original JavaScript implementation of the c-factor BRDF algorithm. The GEE script used is available at this [link](https://code.earthengine.google.com/dadb5bffce05dc52f282b9e0688acd79).

Second, the same scene was fed through the Python implementation using the demonstration script `examples/usage_example.py`. 

The original image and the two BRDF-adjusted images can be visualized in the GEE playground using this [link](https://code.earthengine.google.com/1752f23fda8e3ebbeff5624c618761cb). If you use the inspector tool, you can see that at each pixel, the values of each band are the same for both images. Similarly, the pixel-level percent difference is 0 everywhere... which is a good sign!

## Installation

### A. Simple Option

Simply download or copy the `cfactor.py` script and place it wherever you need it in your local code base or project.

### B. Install from GitHub
1. Clone or download this repository.
2. (Optional) Install it as a Python package:
   ```bash
   pip install .
3. Then Import the module in your scripts and use the `apply_cfactor_brdf_correction` function as needed.