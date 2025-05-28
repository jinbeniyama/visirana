#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Show QE (quantum efficiency).
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  
import astropy.io.fits as fits


def extract_qe_template(fi):
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

    # 0: wavelength in m
    # 1: efficiency
    data = hdu1.data

    arr = np.array(data, dtype=[("w_m", "f8"), ("eff", "f8")])
    df = pd.DataFrame(arr)
    return df


if __name__ == "__main__":
    parser = ap(description="Plot quantum efficiency.")
    parser.add_argument(
        "fi", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "--out", type=str, default="QE.jpg",
        help="Output filename")
    args = parser.parse_args()
    
    # Extract
    df = extract_qe_template(args.fi)
    
    ## Convert [m] to [micron]
    df["w"] = df["w_m"]*1e6
    
    # Plot
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Efficiency")
    ax.set_xlim([8, 14])
    ax.set_ylim([0, 1])
    ax.plot(df["w"], df["eff"], label="Quantum efficiency", color="black")
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
