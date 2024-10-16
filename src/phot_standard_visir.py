#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Do aperture photometry for stacked fits of standard stars taken with VLT/VISIR.
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits
import pandas as pd 
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import sep
from calcerror import adderr

from visirana.util import cut4VISIR, Perr, get_VISIR_standard
from visirana.visualization import plot_growth


    
def phot_ann_stan(
    f1_list, f2_list, pos_list, sign_list, rad, rnin, rout, gain):
    """
    Poisson error is not included in the reported errors.

    Parameters
    ----------
    Return
    ------
    """
    # TODO: 
    # Median stacking?

    N1 = len(f1_list)
    N2 = len(f2_list)
    assert N1==N2, "Check arguments."
    Npos = len(pos_list)
    assert Npos == 4, "Check arguments."

    f1stack_list, f2stack_list, substack_list = [], [], []
    for f1, f2 in zip(f1_list, f2_list):
        hdul_1 = fits.open(f1)
        hdu3_1 = hdul_1[3]
        image_1 = hdu3_1.data
        image_1 = cut4VISIR(image_1)

        #print(f"  Median of Image 1 (original)  : {np.median(image_1):.2f}")
        #image_1 -= np.median(image_1)
        #print(f"  Median of Image 1 (subtracted): {np.median(image_1):.2f}")
        
        hdul_2 = fits.open(f2)
        hdu3_2 = hdul_2[3]
        image_2 = hdu3_2.data
        image_2 = cut4VISIR(image_2)
        #print(f"  Median of Image 2 (original)  : {np.median(image_2):.2f}")
        #image_2 -= np.median(image_2)
        #print(f"  Median of Image 2 (subtracted): {np.median(image_2):.2f}")

        # Subtract
        image_sub = image_1 - image_2
        #print(f"  Median of Image sub (original)  : {np.median(image_sub):.2f}")
        #image_sub -= np.median(image_sub)
        #print(f"  Median of Image sub (subtracted): {np.median(image_sub):.2f}")

        print("")
        # Save 
        f1stack_list.append(image_1)
        f2stack_list.append(image_2)
        substack_list.append(image_sub)

    #imagesub = np.median(substack_list, axis=0) 
    # Stacking 
    imagesub = np.sum(substack_list, axis=0)

    stdsub = np.std(imagesub)

    
    f_list, ferr_list = [], []
    for idx, (xc,yc) in enumerate(pos_list):
        print(f"  Do photometry at ({xc}, {yc}) with rad={rad}")

        sign = sign_list[idx]
        # Cut image to avoid negative values
        if sign == -1:
            img = -1 * imagesub
        else:
            img = imagesub

        # Photometry =======================================  
        # 1. Estimate of background 'noise' using global image:
        # TODO: Update estimate of noise
        bkg = sep.Background(img)
        bgerr_pix = bkg.globalrms

        # 2. Barycenter
        # sigma*4 = 20 pix
        sigma = 5
        xwin, ywin, wflag = sep.winpos(img, xc, yc, sigma)

        # 2. Background error in aperture with local background subtraction with gain = None
        # TODO:
        #    i.e., reported error is simply (apreture area)**2 x bkgerr_pix
        #    Note: If you use sum_circle w/ bkgann, the reported error is 
        #          sum of (apreture area)**2 x bkgerr_pix and some error 
        #          even when gain is None......?
        _, fluxerr, _ = sep.sum_circle(
            img, [xwin], [ywin], r=rad, err=bgerr_pix, gain=None)

        # 3. Aperture flux with local background subtraction with gain=None
        flux, _, _ = sep.sum_circle(
            img, [xwin], [ywin], r=rad, bkgann=(rin, rout), err=bgerr_pix, gain=None)
        flux, fluxerr = float(flux), float(fluxerr)
        # 4. Do not add poisson error to measure the sensitivity
        print(f"    Flux {flux:.2f}+-{fluxerr:.2f}")
        f_list.append(flux)
        ferr_list.append(fluxerr)
        # Photometry =======================================  
        

    # Average
    f = np.sum(f_list)
    ferr = adderr(ferr_list)
    SNR = f/ferr
    print(f"  Summed flux {f:.2f}+-{ferr:.2f} (S/N={SNR:.1f})")
    return f, ferr


def ana_standard(
    f1, f2, gain, pos_list, sign_list, rad_list, annin, annout, df_stan, SNR_lim, ttotal_lim):

    # For standard
    Nnod = 2
    # Extract exposure time of the a fits
    hdu_hdr = fits.open(f1)[0].header
    # This is not an effective exposure time, but a sequence time
    texp_fits = hdu_hdr["HIERARCH ESO DET SEQ1 EXPTIME"] 
    hdu_hdr3 = fits.open(f1)[3].header
    # Number of exposure in 1 chopping
    Nexp_nod = hdu_hdr3["HIERARCH ESO DET FRAM NDIT"]
    # Exposure time per integration
    texp1 = hdu_hdr3["EXPTIME"]
    # This is total exposure time for standard
    ttotal_standard = Nexp_nod*texp1*Nnod

    # TODO: Check
    key_ttotal = "HIERARCH ESO SEQ TIME "
    ttotal_standard = hdu_hdr[key_ttotal]

    # Keyword to extract object like HD145897
    key_obj = "HIERARCH ESO OBS TARG NAME"
    obj = hdu_hdr[key_obj]

    # Keyword to extract filter
    key_fltr = "HIERARCH ESO INS FILT1 NAME"
    fltr = hdu_hdr[key_fltr]
    df_s = df_stan[df_stan["obj"]==obj]
    df_s = df_s.reset_index(drop=False)
    # In Jy
    flux_standard = df_s[fltr][0]

    # Photometry
    #   Pixel scale: 0.045 arcsec/pix 
    #   8 arcsec = 177.77 pix
    #   12 arcsec = 266.667
    #   20 pix = 0.9 arcsec

    flux_list, fluxerr_list, SNR_list, flux_lim_list = [], [], [], []
    for radius in rad_list:
        flux, fluxerr  = phot_ann_stan(
            [f1], [f2], pos_list, sign_list, radius, annin, annout, gain)
        SNR_standard = flux/fluxerr
        
        # Sensitivity estimate
        # SNR = C * t**0.5 * flux
        # Determine C from the flux of Standard star
        C = SNR_standard/flux_standard/(ttotal_standard)**0.5
        # Derive a upper limit of a flux to reach SNR_lim
        flux_lim = SNR_lim/C/(ttotal_lim)**0.5
        flux_lim_mJy = 1000*flux_lim
        print(f"  Sensitivity limit (SNR={SNR_lim}) {flux_lim_mJy:.2f} mJy")

        flux_list.append(flux)
        fluxerr_list.append(fluxerr)
        SNR_list.append(SNR_standard)
        # In mJy
        flux_lim_list.append(flux_lim_mJy)

    # Make a dataframe
    df = pd.DataFrame({
        "flux":flux_list,
        "fluxerr":fluxerr_list,
        "snr":SNR_list,
        "fluxlim":flux_lim_list,
        "radius":rad_list,
        })
    df["band"] = fltr
    df["fluxcat"] = flux_standard
    df["texp"] = ttotal_standard
    return df



    

if __name__ == "__main__":
    parser = ap(description="Make a stacked image for VLT/VISIR.")
    parser.add_argument(
        "f1", type=str,
        help="Fits for nodding 1")
    parser.add_argument(
        "f2", type=str,
        help="Fits for nodding 2")
    parser.add_argument(
        "--pos", type=int, nargs="*",
        help="Coordinates")
    parser.add_argument(
        "--sign", type=int, nargs="*",
        help="Signs of signals")
    parser.add_argument(
        "--rad0", type=int, default=1,
        help="Minimum aperture radius in pix.")
    parser.add_argument(
        "--rad1", type=int, default=20,
        help="Maximum aperture radius in pix.")
    parser.add_argument(
        "--rin", type=int, default=20,
        help="Inner annulus in pix.")
    parser.add_argument(
        "--rout", type=int, default=30,
        help="Outer annulus in pix.")
    parser.add_argument(
        "--nodtype", default="pallarel",
        help="Nodding type (parallel or perpendicular)")
    parser.add_argument(
        "--out_photres", type=str, default="photres.txt",
        help="output file")
    parser.add_argument(
        "--tlim", type=float, default=3600., 
        help="Effective exposure time in s to calculate limiting flux")
    parser.add_argument(
        "--snrlim", type=float, default=1.,
        help="SNR of limiting flux")
    parser.add_argument(
        "--sigma", type=int, default=3,
        help="Dynamic range of the image")
    parser.add_argument(
        "--out", type=str, default=None,
        help="output image")
    args = parser.parse_args()
    
    # nodding 1
    f1 = args.f1
    # nodding 2
    f2 = args.f2
    #sub_nodpair([f1], [f2], sigma=3)

    ## This is not used in this analysis. 
    ## See VISIR manual
    ## From fits header 20.00 / [e-/ADU] Gain
    ## This is inverse gain. This should be ok as an input of sep.sum_circle.
    #gain = 20
    gain = None


    pos_list = []
    N_pos = int(len(args.pos)/2)
    for n in range(N_pos):
        pos_list.append((args.pos[2*n:2*n+2]))
    sign_list = args.sign

    # Setting similar to VISIR manual
    ## See MANUAL
    # rin, rout = 20, 30
    rin, rout = args.rin, args.rout
    rad_list = np.arange(args.rad0, args.rad1+1, 1)

    # Get VISIR standard star
    df_stan = get_VISIR_standard()

    # Effective integration time of the fits of interest
    tlim = args.tlim
    # Plot limiting flux of SNR of snrlim
    snrlim = args.snrlim
    df_phot = ana_standard(
        f1, f2, gain, pos_list, sign_list, rad_list, rin, rout, df_stan, snrlim, tlim)

    # Output best values
    idx_max = df_phot["snr"].idxmax()
    snr_max = df_phot["snr"][idx_max]
    radius_max = df_phot["radius"][idx_max]
    fluxlim_min = df_phot["fluxlim"][idx_max]
    flux_best = df_phot["flux"][idx_max]
    flux_cat = df_phot["fluxcat"][idx_max]
    texp_standard = df_phot["texp"][idx_max]
    band = df_phot["band"][idx_max]
    print("")
    print("Best values:")
    print(f"    Band         = {band}")
    print(f"    SNR          = {snr_max:.1f}")
    print(f"    Radius       = {radius_max:.1f} pix")
    print(f"    limflux      = {fluxlim_min:.2f} mJy (S/N={snrlim} {tlim} s)")
    print(f"    Cat flux  C1 = {flux_cat} mJy")
    print(f"    Best flux C2 = {flux_best:.1f} ADU")
    print(f"    Exp time  C3 = {texp_standard} s")
    coeff = flux_cat/flux_best*texp_standard
    print(f"    Coeff C1/C2*C3   = {coeff:.4f} Jy s/ADU")

    
    # Plot growth curve
    if args.out:
        plot_growth(df_phot, tlim, snrlim, args.out)
