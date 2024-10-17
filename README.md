# VISIR ANALYSIS  (visirana)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jinbeniyama@oca.eu)

## Overview
Analysys of VLT/VISIR photometry of moving objects could be done in this repository.
The icon is from https://www.eso.org/sci/facilities/paranal/instruments/visir.html.
The official manual can be downloaded from the following link.
https://ftp.eso.org/pub/dfs/pipelines/instruments/visir/visir-pipeline-manual-1.11.pdf


## Installing
```
git clone https://jin_beniyama@bitbucket.org/jin_beniyama/visirana.git
```

## Usage: imaging

```
# Do stacking: make positive and negative images
## A-B-B-A nodding is assumed
nodstack.py --f_list (fits_nodA_1) (fits_nodB_1) (fits_nodB_2) (fits_nodA_2) --out (reduced.fits)
## Nodding pairs can specified.
nodstack.py --fA_list (fits_nodA_1) (fits_nodA_2) --fB_list (fits_nodB_1) (fits_nodB_2) --out (reduced.fits)

# Aperture photometry with local background subtraction using an annulus
## Just output the result and plot them
appphot_visir.py (reduced.fits) --xy_posi 144 126 --xy_nega 144 189 144 63 --rad 7 --rin 9 --rout 11 --out_image (out_image.jpg)
## Save in a output file
appphot_visir.py (reduced.fits) --xy_posi 144 126 --xy_nega 144 189 144 63 --rad 7 --rin 9 --rout 11 --out_image (out_image.jpg) --out_res (out_res.txt) --obj (star1) --band (band1)
```


## Dependencies

This library is depending on `NumPy`, `SciPy`, `SEP`, `Astropy` 
and `Astroquery`.
Scripts are developed on `Python 3.7.10`, `NumPy 1.19.2`, `SciPy 1.6.1`,
`SEP 1.0.3`, `Astropy 4.2` and `Astroquery 0.4.1`.

## LICENCE

This software is released under the MIT License.
