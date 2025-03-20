#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Obtain upper limit flux density.
"""
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import scipy.stats as stats


def check_random(f, coeff, texp):
    """
    Calculate flux density with conversion factor.

    Parameters
    ----------
    f : str
        filename
    coeff : float
        conversion factor in Jy/(ADU s^-1)
    texp : float
        exposure time in s

    Return
    ------
    df : pandas.DataFrame
        dataframe with flux density
    """
    # Conversion factor from count to flux density, coeff [Jy/(ADU s^-1)]
    # i.e.,
    # F_target = coeff*C_target

    df = pd.read_csv(f, sep=" ")
    # Flux [Jy]  = coeff [Jy/(ADU s^-1)]  * Count [AUD] / time
    df["fluxdensity"] = coeff * (df["flux"] / texp)
    # mJy
    df["fluxdensity"] *= 1000
    return df


if __name__ == "__main__":
    parser = ap(description="Determine upper limit of flux density.")
    parser.add_argument(
        "res", type=str,
        help="Photometric result.")
    parser.add_argument(
        "--teff", type=float, 
        help="On source exposure time")
    parser.add_argument(
        "--coeff", type=float, 
        help="Coefficient Jy/(ADU s^-1), output of phot_standard_visir.py")
    parser.add_argument(
        "--nbin", type=int, default=30, 
        help="Number of bins of histogram")
    parser.add_argument(
        "--out", type=str, default=None,
        help="output image")
    args = parser.parse_args()
    
    
    # Coefficient mJy s/ADU
    coeff = args.coeff
    # On source exposure time
    teff = args.teff
    df = check_random(args.res, coeff, teff)
    N0 = len(df)

    # TODO: Check
    # Only positive?
    #df = df[df["flux"] > 0]
    Npos = len(df)
    print(f"N all           = {N0}")
    print(f"N positive flux = {Npos}")
    
    f_max = np.max(df["fluxdensity"])
    print(f"  Maximum flux density {f_max:.1f} mJy")

    # Median, 1, 2, 3 sigma
    # TODO Check
    percentage = [67, 90, 99.7]
    for p in percentage:
        fluxdensity_sigma = stats.scoreatpercentile(df["fluxdensity"], p)
        print(f"  Flux density < {p}% = {fluxdensity_sigma:.4f} mJy")

    # Make histogram
    if args.out:
        Nbin = args.nbin

        fig = plt.figure(figsize=(8, 6))
        
        #ax2 = fig.add_axes([0.15, 0.15, 0.30, 0.8])
        #ax2.set_xlabel("Flux [ADU]")
        #ax2.set_ylabel("N")
        #ax2.hist(df["flux"], bins=Nbin, histtype="step")
        
        ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
        ax.set_xlabel("Flux density [mJy]")
        ax.set_ylabel("N")
        label = f"N={len(df)}\n Maximum: {f_max:.2f} mJy"
        col = "black"
        ax.hist(
            df["fluxdensity"], bins=Nbin, color=col, histtype="step", label=label)
        xmin, xmax = ax.get_xlim()
        vmax = np.max([abs(xmin), abs(xmax)])
        ax.set_xlim([-vmax, vmax])
        #ax.vlines(
        #    fluxdensity_sigma, ymin, ymax, color="red", ls="dashed", 
        #    label=f"3-sigma {fluxdensity_sigma:.1f} mJy")
        ax.legend()
        
        plt.savefig(args.out)
