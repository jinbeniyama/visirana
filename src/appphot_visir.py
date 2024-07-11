#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Do aperture photometry for stacked fits taken with VLT/VISIR.
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits
import pandas as pd 
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
import sep
from myplot import mycolor, mymark
from calcerror import adderr


def main(args):
    """This is the main function called by the `appphot_visir.py` script.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments passed from the command-line as defined below.
    """

    #args.xy
    #args.nodtype
    #args.out
    #args.sigma
    
    # Useless?
    # Sometimes some signals are out of field of view.
    # 4 sources for perpendicular
    #if len(args.xy)==8:
    #    assert args.nodtype=="parallel", "Invalid source coordinates and/or nodding type."
    #elif len(args.xy)==6:
    #    assert args.nodtype=="perpendicular", "Invalid source coordinates and/or nodding type."
    
    # Positive signals
    xposi_list = args.xy_posi[::2]
    yposi_list = args.xy_posi[1::2]
    # Negative signals
    xnega_list = args.xy_nega[::2]
    ynega_list = args.xy_nega[1::2]
    rad, rin, rout = args.rad, args.rin, args.rout
    
    hdu = fits.open(args.fi)
    img = hdu[0].data
    hdr = hdu[0].header
    img = img.byteswap().newbyteorder()
    sigma = args.sigma
    
    kwd_gain = "HIERARCH ESO DET CHIP GAIN"
    if kwd_gain in hdr:
        gain = hdr[kwd_gain]
    else:
        assert False, "Gain is unknown."

    # What is the gain for old fits? such as 2009?
    # gain = 1 (old fits, such as 2011)
    # gain = 20 (new fits, such as 2024)
    print(f"Inverse Gain = {gain} e/ADU")

    # Do aperture photometry ==================================================
    f_list, ferr_list = [], []
    # Positive
    for xc, yc in zip(xposi_list, yposi_list):
        flux, fluxerr, eflag = sep.sum_circle(
            img, [xc], [yc], r=rad, gain=gain, bkgann=(rin, rout))
        flux, fluxerr = float(flux), float(fluxerr)
        print(f"    -> flux = {flux:.2f}+-{fluxerr:.2f}")
        print(f"       S/N={flux/fluxerr:.1f})")
        f_list.append(flux)
        ferr_list.append(fluxerr)

    # Negative
    img = -1 * img
    for xc, yc in zip(xnega_list, ynega_list):
        flux, fluxerr, eflag = sep.sum_circle(
            img, [xc], [yc], r=rad, gain=gain, bkgann=(rin, rout))
        flux, fluxerr = float(flux), float(fluxerr)
        print(f"    -> flux = {flux:.2f}+-{fluxerr:.2f}")
        print(f"       S/N={flux/fluxerr:.1f})")
        f_list.append(flux)
        ferr_list.append(fluxerr)
    img = -1 * img

    # Average
    f_ave = np.sum(f_list)
    ferr_ave = adderr(ferr_list)
    SNR_ave = f_ave/ferr_ave
    print(f"  Averaged flux: {f_ave:.2f}+-{ferr_ave:.2f} (S/N={SNR_ave:.1f})")
    # Do aperture photometry ==================================================

    # Plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    ax.set_xlabel("x [pix]")
    ax.set_ylabel("y [pix]")
    stdsub = np.std(img)
    vmin, vmax = -sigma*stdsub, sigma*stdsub

    # Add apertures
    color_1 = "black"
    color_2 = "white"
    ax.imshow(img, cmap='inferno', vmin=vmin, vmax=vmax)

    for xc, yc in zip(xposi_list, yposi_list):
        # radius
        ax.add_collection(PatchCollection(
            [Circle((xc, yc), rad)], color=color_1, ls="solid",
            lw=2, facecolor="None", label=None)
            )
        # annulus
        ax.add_collection(PatchCollection(
            [Circle((xc, yc), rin)], color=color_1, ls="dashed",
            lw=2, facecolor="None", label=None)
            )
        ax.add_collection(PatchCollection(
            [Circle((xc, yc), rout)], color=color_1, ls="dotted",
            lw=2, facecolor="None", label=None)
            )

    for xc, yc in zip(xnega_list, ynega_list):
        # radius
        ax.add_collection(PatchCollection(
            [Circle((xc, yc), rad)], color=color_2, ls="solid",
            lw=2, facecolor="None", label=None)
            )
        # annulus
        ax.add_collection(PatchCollection(
            [Circle((xc, yc), rin)], color=color_2, ls="dashed",
            lw=2, facecolor="None", label=None)
            )
        ax.add_collection(PatchCollection(
            [Circle((xc, yc), rout)], color=color_2, ls="dotted",
            lw=2, facecolor="None", label=None)
            )

    plt.savefig(args.out_image)



    

if __name__ == "__main__":
    parser = ap(description="Make a stacked image for VLT/VISIR.")
    parser.add_argument(
        "fi", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "--xy_posi", type=int, nargs="*",
        help="List of coordinates of 'positive' signals.")
    parser.add_argument(
        "--xy_nega", type=int, nargs="*",
        help="List of coordinates of 'negative' signals.")
    parser.add_argument(
        "--rad", type=int, default=5,
        help="Aperture radius in pix.")
    parser.add_argument(
        "--rin", type=int, default=7,
        help="Inner annulus in pix.")
    parser.add_argument(
        "--rout", type=int, default=9,
        help="Outer annulus in pix.")
    parser.add_argument(
        "--nodtype", default="pallarel",
        help="Nodding type (parallel or perpendicular)")
    parser.add_argument(
        "--out_photres", type=str, default="photres.txt",
        help="output file")
    parser.add_argument(
        "--sigma", type=int, default=3,
        help="Dynamic range of the image")
    parser.add_argument(
        "--out_image", type=str, default="photres.jpg",
        help="output image")
    args = parser.parse_args()
    
    main(args)
