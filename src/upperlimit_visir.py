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
    # Conversion factor from count to flux density
    # C_target = A*F_target
    #texp_target = 1017
    # A = C_standard * texp_standard/texp_target / F_standard
    # HD145897 on 2024-05-26, B10.7
    # flux density of the star in mJy
    #F_standard  = 6.0774e3
    # exposure time of the star in s
    #texp_standard = 120
    #C_standard = 65773.93
    #df["fluxdensity"] = df["flux"] / C_standard / texp_target * texp_standard * F_standard

    df = pd.read_csv(f, sep=" ")

    df["fluxdensity"] = df["flux"] / texp * coeff
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
        help="Coefficient mJy s/ADU, output of phot_standard_visir.py")
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
    df = df[df["flux"] > 0]
    Npos = len(df)
    print(f"N all           = {N0}")
    print(f"N positive flux = {Npos}")
    
    f_max = np.max(df["fluxdensity"])
    print(f"  Maximum flux density {f_max:.1f} mJy")

    # Median, 1, 2, 3 sigma
    percentage = [67, 90, 99.7]
    for p in percentage:
        fluxdensity_sigma = stats.scoreatpercentile(df["fluxdensity"], p)
        print(f"  Flux density < {p}% = {fluxdensity_sigma:.4f} mJy")

    # Make histogram
    if args.out:
        Nbin = 30

        fig = plt.figure(figsize=(12, 6))
        
        ax2 = fig.add_axes([0.15, 0.15, 0.30, 0.8])
        ax2.set_xlabel("flux [ADU]")
        ax2.set_ylabel("N")
        ax2.hist(df["flux"], bins=Nbin, histtype="step")
        
        ax3 = fig.add_axes([0.65, 0.15, 0.30, 0.8])
        ax3.set_xlabel("flux density [mJy]")
        ax3.set_ylabel("N")
        ax3.hist(df["fluxdensity"], bins=Nbin, histtype="step", label=label)
        ymin3, ymax3 = ax3.get_ylim()
        ax3.vlines(fluxdensity_sigma, ymin3, ymax3, color="red", ls="dashed", 
            label=f"3-sigma {fluxdensity_sigma:.1f} mJy")
        ax3.set_ylim(ymin3, ymax3)
        ax2.legend()
        ax3.legend()
        
        plt.savefig(args.out)
