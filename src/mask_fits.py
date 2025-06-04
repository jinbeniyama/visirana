#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Mask some region or cut some region in a fits file and save it.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  
import astropy.io.fits as fits


if __name__ == "__main__":
    parser = ap(description="Mask or cut fits.")
    parser.add_argument(
        "fi", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "mode", type=str,
        help="mask or cut")
    parser.add_argument(
        "xr", nargs=2, type=int,
        help="X range to be masked.")
    parser.add_argument(
        "yr", nargs=2, type=int,
        help="Y range to be masked.")
    parser.add_argument(
        "--out", type=str, default="masked.fits",
        help="Output filename")
    args = parser.parse_args()
    

    # fits source  (HDU0, header + image data)
    hdu = fits.open(args.fi)
    # Update header
    hdu[0].header.add_history(f"[mask_fits.py] Mask added.")

    # Mask region of interest
    xmin, xmax = args.xr
    ymin, ymax = args.yr

    if args.mode == "mask":
        hdu[1].data[ymin:ymax, xmin:xmax] = 0
        hdu[2].data[ymin:ymax, xmin:xmax] = 0
        hdu[3].data[ymin:ymax, xmin:xmax] = 0
   
    elif args.mode == "cut":
        hdu[1].data = hdu[1].data[ymin:ymax, xmin:xmax]
        hdu[2].data = hdu[2].data[ymin:ymax, xmin:xmax]
        hdu[3].data = hdu[3].data[ymin:ymax, xmin:xmax]

    out = args.out
    hdu.writeto(out, overwrite=True)
