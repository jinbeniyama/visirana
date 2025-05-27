# VISIR ANALYSIS  (visirana)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jinbeniyama@oca.eu)

## Overview
Analysys of VLT/VISIR photometry and spectroscopy could be done in this repository.
The official manual can be downloaded from [here](https://ftp.eso.org/pub/dfs/pipelines/instruments/visir/visir-pipeline-manual-1.11.pdf).


## Installing
```
git clone https://jin_beniyama@bitbucket.org/jin_beniyama/visirana.git
```

## Usage

### Stacking (for both imaging and spectroscopy)
```
# Do stacking: make positive and negative images
## A-B-B-A nodding is assumed
nodstack.py --f_list (fits_nodA_1) (fits_nodB_1) (fits_nodB_2) (fits_nodA_2) --out (reduced.fits)
## Nodding pairs can specified.
nodstack.py --fA_list (fits_nodA_1) (fits_nodA_2) --fB_list (fits_nodB_1) (fits_nodB_2) --out (reduced.fits)
```


### Aperture photometry (for imaging)
```
# Aperture photometry with local background subtraction using an annulus
## Just output the result and plot them
phot_visir.py (reduced.fits) --xy_posi 144 126 --xy_nega 144 189 144 63 --rad 7 --rin 9 --rout 11 --out_image (out_image.jpg)
## Save in a output file
phot_visir.py (reduced.fits) --xy_posi 144 126 --xy_nega 144 189 144 63 --rad 7 --rin 9 --rout 11 --out_image (out_image.jpg) --out_res (out_res.txt) --obj (star1) --band (band1)

# Calculate noise level from aperture photometry of standard star
## Only for perpendicular nodding 1 and 2
phot_standard_visir.py (fits_nod1) (fits_nod2) --tlim (exposure time to calculcate sensitivity) --pos  x1 y1 x2 y2 x3 y3 x4 y4 --sign -1 1 1 -1 --out (outputfile) --rad0 (inner radius in pix) --rad1 (maximum radius in pix)

# Do photometry using random apertures 
## Background noise is not saved. Annuli are used only to estimate background level
phot_random_visir.py (fits) -N (Number of apertures) --rad (aperture radius) --rin (inner edge of annulus) --rout (outer edge of annulus) --out_photres (photometry result)
```

### Extract 1-D spectrum (for spectroscopy)
This process is interactive. 

``` 
# Extract with argmax
spec_visir.py (fits after stacking) argmax
```

``` 
# Extract with 1st-order moment
spec_visir.py (fits after stacking) moment
```

### Wavelength calibration (for spectroscopy)
In prep.


## Dependencies

This library is depending on `NumPy`, `SciPy`, `SEP`, `Astropy` 
and `Astroquery`.
Scripts are developed on `Python 3.7.10`, `NumPy 1.19.2`, `SciPy 1.6.1`,
`SEP 1.0.3`, `Astropy 4.2` and `Astroquery 0.4.1`.

## LICENCE
This software is released under the MIT License.
