#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Extract fluxes from fits file after reduction using esoreflex pipeline.
"""
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import astropy.io.fits as fits
from astropy.time import Time
import sep


def extract_VISIR_photometry(fi):
    """Perform photometry with calibrated data using pipeline.

    Parameter
    ---------
    fi : str
        fits file contains qe.

    Return
    ------
    df : pandas.DataFrame
        dataframe with wavelength and efficiency
    """
    hdu = fits.open(fi)
    N_hdu = len(hdu)
    hdu1 = hdu[0]
    fltr = hdu[0].header["HIERARCH ESO INS FILT1 NAME"]
    # J12.2 to 12.2
    wave = fltr[1:]
    obj = hdu[0].header["OBJECT"]
    mjd = hdu[0].header["MJD-OBS"]
    jd = Time(str(mjd), format='mjd', scale='utc').jd
    # hdu0: Flux [Jy]
    # hdu1: Flux err [Jy]
    # hdu2: Weight
    # hdu3: Bad pixel

    img = hdu1.data
    img = img.astype(img.dtype.newbyteorder('='))
    ny, nx = img.shape
    xc, yc = nx/2, ny/2
    sigma = 10

    xwin, ywin, wflag = sep.winpos(img, xc, yc, sigma)
    # VLT PIPELINE is using R=30?
    radii = np.arange(5, 41, 1)

    flux, fluxerr, _ = sep.sum_circle(
        img, [xwin], [ywin], r=radii, err=1, gain=None)
    
    print(f"{obj} in {fltr}-band micront at {jd}")
    for f,fe,r in zip(flux, fluxerr, radii):
        print(f"flux = {f:.4f}+-{fe:.4f} (rad={r})")
    print("")

    # Flux in Jy
    df = pd.DataFrame(dict(flux=flux, fluxerr=fluxerr, rad=radii))
    # TODO: Check
    # Assume 10 % error
    df["fluxerr"] = df["flux"]*0.1
    # Wavelength in micron
    df["wavelength"] = float(wave)
    # TODO: Check
    # Consider light-time correction, exposure time, etc.
    df["jd"] = jd
    # TODO: Consider
    df["cflag"] = 0
    df["memo"] = f"VLT_VISIR_{obj}"
    df["code"] = 309
    return df


if __name__ == "__main__":
    parser = ap(description="Extract flux.")
    parser.add_argument(
        "fi", type=str, nargs="*",
        help="Fits file after reduction with esoreflex pipeline.")
    parser.add_argument(
        "--rad", type=int, default=None,
        help="Aperture radius")
    parser.add_argument(
        "--out", type=str, default="flux.txt",
        help="Output filename")
    args = parser.parse_args()
    
    df_list = []
    # Extract
    for fi in args.fi:
        df = extract_VISIR_photometry(fi)
        df_list.append(df)
    df = pd.concat(df_list)
    
    if args.rad:
        print(f"Save only rad=={args.rad} for TPM")
        # Extract one radius for TPM
        df = df[df["rad"] == args.rad]
        df = df.reset_index(drop=True)
        col_use = ["jd", "wavelength", "flux", "fluxerr", "code", "cflag", "memo"]
        # Do not include rad
        df = df[col_use]

    df.to_csv(args.out, sep=" ", index=False)
