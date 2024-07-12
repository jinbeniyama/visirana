# VISIR ANALYSIS  (visirana)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[developer mail](mailto:jinbeniyama@oca.eu)

## Overview

Analysys of VLT/VISIR photometry of moving objects could be done in this repository.
The icon is from https://www.eso.org/sci/facilities/paranal/instruments/visir.html.


## Installing
```
git clone https://jin_beniyama@bitbucket.org/jin_beniyama/visirana.git
```

## Usage: imaging

```
[usage]
# Create database (once)
psdb.py dbstart --first
# Create new tables
psdb.py create --table (table name)
# Insert stars to database
psdb.py insert --tablename (table name) --ra (ra in degree) --dec (dec in degree)
--radius (fov radius in degree) --magmin (minimum magnitude) --magmax (maximum magnitude)
# Check stars in database
psdb.py extract --tablename (table name) --ra (ra in degree) --dec (dec in degree)
--radius (fov radius in degree) --magmin (minimum magnitude) --magmax (maximum magnitude)

[example]
# Create 2021DX1 table.
# Observed locations are below.
psdb.py create --table _2021DX1
psdb.py insert --table _2021DX1 --ra 208.87 --dec 44.68 --radius 0.2 --magmin 12 --magmax 21
```


## Dependencies

This library is depending on `NumPy`, `SciPy`, `SEP`, `Astropy` 
and `Astroquery`.
Scripts are developed on `Python 3.7.10`, `NumPy 1.19.2`, `SciPy 1.6.1`,
`SEP 1.0.3`, `Astropy 4.2` and `Astroquery 0.4.1`.

## LICENCE

This software is released under the MIT License.
