#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Extract 1-d spectrum from stacked fits taken with VLT/VISIR.

Excellent Reference
-------------------
https://learn.astropy.org/tutorials/1-SpectroscopicTraceTutorial.html
"""
from argparse import ArgumentParser as ap
import os
import sys
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from astropy.modeling.polynomial import Polynomial1D
from astropy.modeling.fitting import LinearLSQFitter
from myplot import mycolor


def extract_spec_argmax(img):
    """Extract bright points with argmax.

    Parameter
    ---------
    img : 2d array-like
        2d image of interest

    Returns
    -------
    xvals, yvals : array-like
        x values and extracted argmax values
    """
    ny, nx= img.shape
    yvals = np.argmax(img, axis=0)
    xvals = np.arange(nx)
    return xvals, yvals


def extract_spec_moment(img):
    """Extract bright points with 1st-order moment.

    Parameter
    ---------
    img : 2d array-like
        2d image of interest

    Returns
    -------
    xvals, yvals : array-like
        x values and extracted intensity-weighted mean positions
    """
    ny, nx= img.shape
    ## Try to find the spine to trace using weighted
    # we use a cutout around the traced line, so the Y-values are from that cutout
    # the `repeat` command here is used to extend our Y-axis position values, which are 425, 426, ... 475
    # along the X-direction.  The indexing with [:, None] adds a "dummy" axis along the second (x) dimension,
    # then `repeat` copies our Y-axis values.  The resulting array has the same shape as our weight array,
    # which is image_array[425:475, :] minus the median

    ymin, ymax = 0, ny
    xvals = np.arange(nx)
    yaxis = np.repeat(np.arange(ymin, ymax)[:,None], img.shape[1], axis=1)
    background = np.median(img)
    # Background is already subtracted in most cases......?
    # background = 0
    print(f"  Median background subtracted {background:.2f} ADU")
    # moment 1 is the data-weighted average of the Y-axis coordinates
    yvals = np.average(yaxis, axis=0, weights=img[ymin:ymax,:] - background)
    return xvals, yvals


def extract_effective_spec(
    fi, mode, xmin=None, xmax=None, ymin=None, ymax=None, fit_degree=2):
    """Extract effective 1d spectrum.

    Parameters
    ----------
    fi : str
        preprocessed fits file (i.e., hdu[0].data is stacked image)
    mode : str
        "argmax" (arguments of the maxima) or "moment" (1st-order moment)
    xmin, xmax : int, optional
        minimum and maximum x values of interest
    ymin, ymax : int, optional
        minimum and maximum y values of interest
    fit_degree : int, optional
        spectrum is fitted [fit_degree]th-order polynomial 

    Returns
    -------
    xeff, yeff : array-like
        x and y values
    yfit : array-like
        fitted y values
    """

    # Fits file of N-band spectroscopy of standard star after stacking
    hdu = fits.open(fi)
    img = hdu[0].data
    ny, nx = img.shape

    if xmin is None:
        xmin = 0
    if xmax is None:
        xmax = nx
    if ymin is None:
        ymin = 0
    if ymax is None:
        ymax = ny

    # All
    if mode == "argmax":

        while True:
            print(f"xmin, xmax, ymin, ymax = {xmin}, {xmax}, {ymin}, {ymax}")
            xv, yv = extract_spec_argmax(img)
            # Remove outliers with ymin and ymax
            # Should be iterative
            bad_x = (xv < xmin) | (xv > xmax)
            bad_y = (yv < ymin) | (yv > ymax)
            xv_bad_y, yv_bad_y = xv[bad_y], yv[bad_y]
            xv_bad_x, yv_bad_x = xv[bad_x], yv[bad_x]

            fig = plt.figure(figsize=(20, 4))
            ax1 = fig.add_axes([0.1, 0.1, 0.20, 0.8])
            ax2 = fig.add_axes([0.32, 0.1, 0.20, 0.8])
            ax3 = fig.add_axes([0.54, 0.1, 0.20, 0.8])
            ax4 = fig.add_axes([0.76, 0.1, 0.20, 0.8])

            ax1.scatter(xv, yv, color="gray", marker='x', label=f"All N={len(xv)}")
            ax1.set_ylabel("Argmax trace data")
            ax1.set_xlabel("x [pixel]")

            ax2.scatter(xv, yv, color=mycolor[1], marker="x")
            ax2.scatter(xv[bad_y], yv[bad_y], color="red", marker="x", label=f"Bad y N={len(xv[bad_y])}")
            ax2.set_xlabel("x [pixel]")

            ax3.scatter(xv, yv, color=mycolor[1], marker="x")
            ax3.scatter(xv[bad_x], yv[bad_x], color="red", marker="x", label=f"Bad x N={len(xv[bad_x])}")
            ax3.set_xlabel("x [pixel]")

            # Effective data
            xeff, yeff = xv[~(bad_x|bad_y)], yv[~(bad_x|bad_y)]

            # Nth-order polynomial
            polymodel = Polynomial1D(degree=fit_degree)
            linfitter = LinearLSQFitter()
            y_fit = linfitter(polymodel, xeff, yeff)
            yfit = y_fit(xeff)
            ax4.scatter(xeff, yeff, color=mycolor[1], s=3, label=f"Selected N={len(xeff)}")
            ax4.scatter(xeff, yfit, color=mycolor[2], marker="s", s=1, label=f"Fit")
            ax4.set_xlabel("x [pixel]")

            for ax in fig.axes:
                ax.legend()
            plt.show(block=False)

            ans = input("Change limits? (y/n): ").strip().lower()
            if ans != 'y':
                break
            try:
                xmin = int(input(f"xmin? (default: {xmin}) ") or xmin)
                xmax = int(input(f"xmax? (default: {xmax}) ") or xmax)
                ymin = int(input(f"ymin? (default: {ymin}) ") or ymin)
                ymax = int(input(f"ymax? (default: {ymax}) ") or ymax)
            except ValueError:
                print("Invalid input. Exiting.")
                break
            plt.close()
        plt.close()

        return xeff, yeff, yfit

    elif mode == "moment":
        while True:
            print(f"xmin, xmax, ymin, ymax = {xmin}, {xmax}, {ymin}, {ymax}")

            # Original
            xv0, yv0 = extract_spec_moment(img)
            img_cut = img[ymin:ymax, xmin:xmax]
            xv, yv = extract_spec_moment(img_cut)

            # Nth-order polynomial
            polymodel = Polynomial1D(degree=fit_degree)
            linfitter = LinearLSQFitter()
            y_fit = linfitter(polymodel, xv, yv)
            yfit = y_fit(xv)

            # Plot
            fig = plt.figure(figsize=(20, 4))
            ax1 = fig.add_axes([0.1, 0.1, 0.35, 0.8])
            ax2 = fig.add_axes([0.55, 0.1, 0.35, 0.8])

            ax1.scatter(xv0, yv0, color="gray", marker='x', label=f"All N={len(xv0)}")
            ax1.set_ylabel("y [pixel]")
            ax1.set_xlabel("x [pixel]")

            ax2.scatter(xv, yv, color=mycolor[1], marker="x", s=3, label=f"Selected N={len(xv)}")
            ax2.scatter(xv, yfit, color=mycolor[2], marker="x", s=1, label="Fit")
            ax2.set_xlabel("x [pixel]")

            for ax in fig.axes:
                ax.legend()

            plt.show(block=False)

            ans = input("Change limits? (y/n): ").strip().lower()
            if ans != 'y':
                break
            try:
                xmin = int(input(f"xmin? (default: {xmin}) ") or xmin)
                xmax = int(input(f"xmax? (default: {xmax}) ") or xmax)
                ymin = int(input(f"ymin? (default: {ymin}) ") or ymin)
                ymax = int(input(f"ymax? (default: {ymax}) ") or ymax)
            except ValueError:
                print("Invalid input. Exiting.")
                break
            plt.close()
        plt.close()
        
        # This is mandatory since the image is cut above
        xeff = [x+xmin for x in xv]
        yeff = [y+ymin for y in yv]
        yfit = [y+ymin for y in yfit]
        return xeff, yeff, yfit



def main(args):
    """This is the main function called by the `spec_visir.py` script.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments passed from the command-line as defined below.
    """

    fi   = args.fi
    mode = args.mode

    # Extract effective spectrum with fitting
    # TODO: Remove args.mode from input arguments and plot both?
    xv, yv, yfit = extract_effective_spec(
        fi, mode, 
        xmin=args.xmin, xmax=args.xmax, 
        ymin=args.ymin, ymax=args.ymax
        )

    # Plot two results on the images
    # Choose spectrum (method)
    # Save spectrum
    # x, ADU

    hdu = fits.open(fi)
    img = hdu[0].data
    ny, nx = img.shape
    
    if args.xr:
        x0, x1 = args.xr
    else:
        x0, x1 = 0, nx
    if args.yr:
        y0, y1 = args.yr
    else:
        y0, y1 = 0, ny
    if args.vr:
        v0, v1 = args.vr
    else:
        v0, v1 = 0, 100
    
    #fig = plt.figure(figsize=(30, 4))
    #ax1 = fig.add_axes([0.05, 0.1, 0.17, 0.8])
    #ax2 = fig.add_axes([0.24, 0.1, 0.17, 0.8])
    #ax3 = fig.add_axes([0.43, 0.1, 0.17, 0.8])
    #ax4 = fig.add_axes([0.62, 0.1, 0.17, 0.8])
    #ax5 = fig.add_axes([0.81, 0.1, 0.17, 0.8])
    #
    #for ax in [ax1, ax2, ax3, ax4, ax5]:
    #    ax.set_xlabel("x [pixel]")
    #ax1.set_ylabel("y [pixel]")
    #
    #for ax in [ax1, ax2]:
    #    ax.imshow(img, vmin=v0, vmax=v1)
    #    ax.scatter(
    #        xv, yfit, color="red", marker="x", s=1, label=f"Fit ({mode})")
    #    ax.legend().get_frame().set_alpha(1.0)
    #
    ## Change aspect
    #ax2.set(aspect=5)
    #ax2.set_xlim([x0, x1])
    #ax2.set_ylim([y0, y1])
    #
    ## Photometry region
    #N_phot = args.N_phot
    #yl = [y - N_phot for y in yfit]
    #yu = [y + N_phot for y in yfit]
    #ax.fill_between(
    #    xv, yl, yu, color='orange', alpha=0.5)
    #
    ## start by taking +/- Npix pixels
    #cutouts = np.array([img[int(yv)-N_phot:int(yv)+N_phot, ii]
    #                    for yv, ii in zip(yfit, xv)])
    #ax3.imshow(cutouts.T)
    #ax3.set_aspect(10)
    #
    #mean_trace_profile = (cutouts).mean(axis=0)
    #trace_profile_xaxis = np.arange(-N_phot, N_phot)
    #ax4.plot(trace_profile_xaxis, mean_trace_profile)
    #ax4.set_xlabel("Distance from center")
    #ax4.set_ylabel("Average source profile")
    #average_spectrum = (cutouts).mean(axis=1)
    #ax5.plot(average_spectrum)
    #plt.show()
    #plt.close()

    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1.2])
    
    # Panel 1: Original image + fit
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.imshow(img, vmin=v0, vmax=v1)
    ax1.scatter(xv, yfit, color="red", marker="x", s=1, label=f"Fit ({mode})")
    ax1.set_xlabel("x [pixel]")
    ax1.set_ylabel("y [pixel]")
    ax1.legend(framealpha=1.0)
    ax1.set_title("Original image with fit")
    
    # Panel 2: Zoomed fit with aspect ratio
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.imshow(img, vmin=v0, vmax=v1)
    ax2.scatter(xv, yfit, color="red", marker="x", s=1, label=f"Fit ({mode})")
    ax2.set(aspect=5)
    ax2.set_xlim([x0, x1])
    ax2.set_ylim([y0, y1])
    yl = [y - args.N_phot for y in yfit]
    yu = [y + args.N_phot for y in yfit]
    ax2.fill_between(xv, yl, yu, color='orange', alpha=0.5)
    ax2.set_xlabel("x [pixel]")
    ax2.legend(framealpha=1.0)
    ax2.set_title("Zoomed region with photometry area")
    
    # Panel 3: Cutouts stack
    ax3 = fig.add_subplot(gs[1, :])
    cutouts = np.array([img[int(yv)-args.N_phot:int(yv)+args.N_phot, ii]
                        for yv, ii in zip(yfit, xv)])
    ax3.imshow(cutouts.T, aspect='auto', cmap='gray')
    ax3.set_xlabel("x [pixel]")
    ax3.set_ylabel("y from the center [pixel]")
    ax3.set_title("After correction")
    
    # Panel 4: Average source profile
    ax4 = fig.add_subplot(gs[2, 0])
    mean_trace_profile = cutouts.mean(axis=0)
    trace_profile_xaxis = np.arange(-args.N_phot, args.N_phot)
    ax4.plot(trace_profile_xaxis, mean_trace_profile)
    ax4.set_xlabel("Distance from center")
    ax4.set_ylabel("Average flux [ADU]")
    ax4.set_title("Average trace profile")
    
    # Panel 5: Spectrum
    ax5 = fig.add_subplot(gs[2, 1])
    average_spectrum = cutouts.mean(axis=1)
    ax5.plot(average_spectrum)
    ax5.set_xlabel("x [pixel]")
    ax5.set_ylabel("Average flux [ADU]")
    ax5.set_title("Extracted spectrum")
    
    plt.tight_layout()
    plt.show(block=False)

    ans = input("Output spectrum? (y/n): ").strip().lower()
    if ans != 'y':
        plt.close()
    try:
        out = args.out_specres
        df = pd.DataFrame(
            {"x":xv, "flux":average_spectrum})
        # Save 
        df.to_csv(out, sep=" ", index=False)

    except ValueError:
        print("Not saved. Exiting.")
    plt.close()


if __name__ == "__main__":
    parser = ap(description="Extract 1-d spectrum from VISIR spectroscopy data.")
    parser.add_argument(
        "fi", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "mode", type=str,
        help="argmax or moment")
    parser.add_argument(
        "--xmin", type=int, default=None,
        help="Minimum x value")
    parser.add_argument(
        "--xmax", type=int, default=None,
        help="Maximum x value")
    parser.add_argument(
        "--ymin", type=int, default=None,
        help="Minimum y value")
    parser.add_argument(
        "--ymax", type=int, default=None,
        help="Maximum y value")
    parser.add_argument(
        "--fit_degree", type=int, default=2,
        help="Degree of polynomial")
    parser.add_argument(
        "--xr", type=int, nargs=2, default=None,
        help="x range of the figure")
    parser.add_argument(
        "--yr", type=int, nargs=2, default=None,
        help="y range of the figure")
    parser.add_argument(
        "--vr", type=int, nargs=2, default=None,
        help="value range of the figure")
    parser.add_argument(
        "--N_phot", type=int, default=10,
        help="Photometry radius in pix.")
    parser.add_argument(
        "--out_specres", type=str, default="specres.txt",
        help="output file")
    args = parser.parse_args()
    
    main(args)
