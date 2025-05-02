#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Show image and save it.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib.pyplot as plt  
import astropy.io.fits as fits

from visirana.visualization import mycolor, myls, mymark
from visirana.util import obtain_winpos, remove_background2d, fit_psf, radial_profile

if __name__ == "__main__":
    parser = ap(description="Plot PSF.")
    parser.add_argument(
        "fi", nargs="*", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "-x", nargs="*", type=float, 
        help="object x in pixel coordinate")
    parser.add_argument(
        "-y", nargs="*", type=float,
        help="object y in pixel coordinate")
    parser.add_argument(
        "--obj", nargs="*", type=str,
        help="object name")
    parser.add_argument(
        "--fittype", type=str, default="Moffat",
        help="Fitting type Gauss or Moffat")
    parser.add_argument(
        "-w", "--width", type=int, default=15,
        help="widht of the cutting region")
    parser.add_argument(
        "-r", "--radius", type=int, default=10,
        help="radius to obtain winpos")
    parser.add_argument(
        "--p_scale", type=float, default=0.35,
        help="pixel scale in arcsec/s")
    parser.add_argument(
        "--norm", action="store_true", default=False,
        help="Normalize brightness by the model peak")
    parser.add_argument(
        "--log", action="store_true", default=False,
        help="Log scale brightness")
    parser.add_argument(
        "--out", type=str, default="psf.png",
        help="Output filename")
    parser.add_argument(
        "--outdir", type=str, default=".",
        help="output directory")
    args = parser.parse_args()
    

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    N_fits = len(args.fi)
    N_obj  = len(args.x)
    
    # Picel scale arcsec/pix
    p_scale = 0.045
    oneside = True
  
    image_list, bgerr_list = [], []
    for infits in args.fi:
        # fits source  (HDU0, header + image data)
        src = fits.open(infits)[0]
        # fits header
        hdr = src.header
        # Gain
        gain = hdr["HIERARCH ESO DET CHIP GAIN"]
        # 2-d image
        image = src.data.byteswap().newbyteorder()
        # BG subtraction
        image, bg_info = remove_background2d(image, None)
        bgerr = bg_info["rms"]
        image_list.append(image)
        bgerr_list.append(bgerr)


    width = args.width
    radius = args.radius
    ny, nx = image.shape
    print(f"Shape of image nx, ny = {nx}, {ny}")
    
    # Calculate fwhm of the fits with clean stars
    #fwhm_mean, fwhm_std = calc_FWHM_fits(infits, radius)
    #fwhm_mean_arcsec = fwhm_mean*p_scale
    #fwhm_std_arcsec = fwhm_std*p_scale
    #print(f"FWHM of image {fwhm_mean:.3f}+-{fwhm_std:.3f} pix")
    #print(f"          <-> {fwhm_mean_arcsec:.3f}+-{fwhm_std_arcsec:.3f} arcsec")

    
    # PSF Plot ================================================================
    # For filename
    tail = ""
    # Plot fitting result
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes([0.15, 0.13, 0.84, 0.77])
    ax.set_xlabel("Radius [arcsec]")
    if args.norm:
        ax.set_ylabel("Normalized Flux")
    else:
        ax.set_ylabel("Count [ADU]")
    for idx_obj, (x, y) in enumerate(zip(args.x, args.y)):
        print(f"\n  Object {idx_obj} at ({x}, {y})")
        obj = args.obj[idx_obj]
        color = mycolor[idx_obj]
        mark = mymark[idx_obj]
        ls = myls[idx_obj+1]
        image = image_list[idx_obj]
        bgerr = bgerr_list[idx_obj]

        # !! Order is important !!
        # At first, serach barycenter
        print(f"111 {nx}")
        xwin, ywin, _ = obtain_winpos(image, [x], [y], [radius], nx, ny)
        xwin, ywin = xwin[0], ywin[0]
        print("Search baricenter")
        print(f"  {x}, {y} -> {xwin}, {ywin}")

        # Second, cut arrays
        # if xwin = 5.5, ywin = 5.5,
        # 0 1 2 3 4 "5" 6 7 8 9 10 in total 11 pix
        # The same for y-direction
        # 5.5 to 5
        idx_x, idx_y = np.floor(xwin), np.floor(ywin)
        xmin, xmax = idx_x - width, idx_x + width 
        ymin, ymax = idx_y - width, idx_y + width 
        print(f"idx_x, idx_y = {idx_x}, {idx_y}")
        print(f"xmin,  xmax = {xmin}, {xmax}")
        print(f"ymin,  ymax = {ymin}, {ymax}")
        # + 1 is necessary
        image_cut = image[int(ymin):int(ymax+1), int(xmin):int(xmax+1)]
        
        # Now, object is in the central pix 
        center = (width, width)
      
        # Need gain for poisson error calculation
        val, valerr = radial_profile(image_cut, center, bgerr, gain)
        x_fit_pix = np.arange(0, len(val), 1)

        print(f"Shape of radial profile: {len(val)}")

        # For plot
        x_model_pix = np.arange(0, len(val), 0.1)
      
        # Fitting in pixel
        if args.fittype == "Gauss":
            params_G = fit_psf(
                x_fit_pix, val, yerr=valerr, fittype="Gauss", oneside=True)
            if oneside:
                A_G, sigma_G = params_G
                mu_G = 0
            else:
                A_G, mu_G, sigma_G = params_G
            sigma_G_arcsec = sigma_G * p_scale
            FWHM_G = 2.35 * sigma_G
            FWHM_G_arcsec = 2.35 * sigma_G * p_scale
            print(f"Params: A_G, mu_G, sigma_G = {A_G:.1f}, {mu_G:.1f}, {sigma_G:.1f}")
            y_model = A_G*np.exp(-0.5*(x_model_pix)**2/sigma_G**2)
            label = (
                f"Gauss FWHM = {FWHM_G_arcsec:.2f} arcsec")

        if args.fittype == "Moffat":
            # Moffat fit
            params_M = fit_psf(
                x_fit_pix, val, yerr=valerr, fittype="Moffat", oneside=True)
            if oneside:
                A_M, sigma_M, beta_M = params_M
            else:
                A_M, mu_M, sigma_M, beta_M = params_M
            # See Tabeshian+2019 3.1
            FWHM_M = 2*sigma_M*np.sqrt(2**(1/beta_M)-1)
            FWHM_M_arcsec = FWHM_M * p_scale
            print(f"Params: A_M, sigma_M, beta_M = {A_M:.1f}, {sigma_M:.1f}, {beta_M:.1f}")
            y_model = A_M*(1+(x_model_pix)**2/sigma_M**2)**(-beta_M)
            label = (
                f"Moffat FWHM = {FWHM_M_arcsec:.3f} arcsec")


        # Plot ================================================================
        if args.norm:
            norm = np.max(y_model)
        else:
            norm = 1.
        y_model = [x/norm for x in y_model]
        val = [x/norm for x in val]
        valerr = [x/norm for x in valerr]

        # Plot in arcsec
        x_fit_arcsec = [x*p_scale for x in x_fit_pix]
        x_model_arcsec = [x*p_scale for x in x_model_pix]

        ax.plot(x_model_arcsec, y_model, label=label, color=color, ls=ls)


        # Plot Data in arcsec
        ax.scatter(
            x_fit_arcsec, val, label=f"{obj} Obs.", 
            facecolor="None", edgecolor=color, marker=mark, s=100)
        ax.errorbar(
            x_fit_arcsec, val, valerr, color=color, 
            lw=2, marker=" ", capsize=0, ls="None")
        # Plot ================================================================

    ax.grid(which="major", axis="both")
    ax.legend(fontsize=20, loc="upper right")
    if args.log:
        ax.set_yscale("log")
        ax.set_ylim([1e-3, 2])

    out = os.path.join(args.outdir, args.out)
    plt.savefig(out)
    # PSF Plot ================================================================
