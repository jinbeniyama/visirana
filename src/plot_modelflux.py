#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Show model flux(es) of standard star(s).
"""
import os 
import sys
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  
import astropy.io.fits as fits


def extract_modelflux(fi, obj):
    """Extract qe (quantum efficiency).

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
    hdu0 = hdu[0]
    hdu1 = hdu[1]

    for n in range(len(hdu1.data)):
        data = hdu1.data[n]
        # 0: Object name
        # 1: RA in deg
        # 2: DEC in deg
        # 3: wavelength [micron]
        # 4: model flux [W/m^2/m]
        obj_data = data[0]
        if obj_data == obj:

            w = data[3].byteswap().newbyteorder()
            flux = data[4].byteswap().newbyteorder()

            df = pd.DataFrame({"w":w, "flux":flux})
            return df
    return None


if __name__ == "__main__":
    parser = ap(description="Plot model flux of standard star.")
    parser.add_argument(
        "fi", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "obj", nargs="*", type=str,
        help="Objects of interest.")
    parser.add_argument(
        "--out", type=str, default="modelflux.jpg",
        help="Output filename")
    args = parser.parse_args()
    

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Flux density [Jy]")
    ax.set_yscale("log")

    for ob in args.obj:
        df = extract_modelflux(args.fi, ob)
        if df is None:
            print(f"{ob} not found in {args.fi}.")
            sys.exit()

        # Original unit is [micron] and [W/m^2/m] (See 7.2 Spectroscopy Standard Stars of visir-pipeline-manual-1.11.pdf)
        ## Convert [micron] to [m]
        df["w_m"] = df["w"]*1e-6
        ## Convert [W/m^2/m] to [Jy]
        # speed of light [m/s]
        c = 2.998e8
        df["flux_jy"] = (df["w_m"]**2/c * df["flux"]) / 1e-26
        ax.plot(df["w"], df["flux_jy"], label=f"Model {ob}")

    ax.legend()
    plt.show(block=False)

    ans = input("Save figure? (y/n): ").strip().lower()
    if ans != 'y':
        plt.close()
    try:
        plt.savefig(args.out)

    except ValueError:
        print("Not saved. Exiting.")
    plt.close()
