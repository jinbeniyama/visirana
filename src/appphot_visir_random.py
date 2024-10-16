#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Do aperture photometry for stacked fits taken with VLT/VISIR.
Just measure the flux in the aperture with local background subtraction.
"""
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd 
import astropy.io.fits as fits
import sep
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
import matplotlib.pyplot as plt


def main(args):
    """
    This is the main function called by the `appphot_visir_random.py` script.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments passed from the command-line as defined below.
    """

    # Open the fits ===========================================================
    # Gain is not used to calculate just "flux"
    hdu = fits.open(args.fi)
    img = hdu[0].data
    hdr = hdu[0].header
    img = img.byteswap().newbyteorder()
    sigma = args.sigma
    rad, rin, rout = args.rad, args.rin, args.rout

    seed = 0
    np.random.seed(seed)
    
    # Do not use edge
    Npix_y, Npix_x = img.shape
    margin = 5*rad
    x_list = np.random.randint(low=margin, high=Npix_x-margin, size=args.N)
    y_list = np.random.randint(low=margin, high=Npix_y-margin, size=args.N)
    # Open the fits ===========================================================
    

    # Do aperture photometry ==================================================
    f_list = []
    # TODO: Is this correct?
    # Flux should be positive
    for xc, yc in zip(x_list, y_list):
        flux, _, _ = sep.sum_circle(img, [xc], [yc], r=rad, bkgann=(rin, rout))
        if flux < 0:
            img_inv = -1 * img
            flux, _, _ = sep.sum_circle(
                img, [xc], [yc], r=rad, bkgann=(rin, rout))
            
        flux = float(flux)
        print(f"    -> flux = {flux:.2f}")
        f_list.append(flux)
    # Do aperture photometry ==================================================


    # Save photometric results ================================================
    df = pd.DataFrame(dict(flux=f_list))
    df["rad"] = rad
    df["rin"] = rin
    df["rout"] = rout
    out = args.out_photres
    df.to_csv(out, index=False, sep=" ")
    # Save photometric results ================================================


    # Plot and save ===========================================================
    if args.out_image:
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

        for xc, yc in zip(x_list, y_list):
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

        plt.savefig(args.out_image)
    # Plot and save ===========================================================


if __name__ == "__main__":
    parser = ap(description="Put random apertures.")
    parser.add_argument(
        "fi", type=str,
        help="Fits after stacking.")
    parser.add_argument(
        "-N", type=int, default=10,
        help="Number of random aperture")
    parser.add_argument(
        "--rad", type=int, default=5,
        help="Aperture radius in pix.")
    parser.add_argument(
        "--rin", type=int, default=20,
        help="Inner annulus in pix.")
    parser.add_argument(
        "--rout", type=int, default=30,
        help="Outer annulus in pix.")
    parser.add_argument(
        "--out_photres", type=str, default="photres.txt",
        help="output file")
    parser.add_argument(
        "--out_image", type=str, default=None,
        help="output image")
    parser.add_argument(
        "--sigma", type=int, default=3,
        help="Dynamic range of the image")
    args = parser.parse_args()
    
    main(args)
