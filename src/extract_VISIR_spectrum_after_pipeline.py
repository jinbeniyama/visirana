#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Extract spectrum from M....fits.

Note: Emissivity spectrum is obtained by dividing observation by polynomial here, 
which is not physical. So output file does not contain emissivity.
Please do TPMs by yourself.

Output file contains followins columns. Thus, this output can be used for Marco's TPM.
- jd      : Julian day 
- flux    : Flux [Jy]
- fluxerr : Flux error [Jy]
- code    : MPC code 
- cflag   : flag for color term correction
- memo    : Memo
"""
import os 
import sys
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  
import astropy.io.fits as fits
from astropy.modeling.polynomial import Polynomial1D
from astropy.modeling.fitting import LinearLSQFitter
from astropy.time import Time

from visirana.util import telluric_features

def extract_1dspec(fi):
    """Extract calibrated spectrum.

    Parameter
    ---------
    fi : str
        fits file contains qe.

    Return
    ------
    df : pandas.DataFrame
        dataframe with wavelength [micron] and efficiency
    """
    hdu = fits.open(fi)
    N_hdu = len(hdu)
    hdu1 = hdu[1]
    obj = hdu[0].header["OBJECT"]
    mjd = hdu[0].header["MJD-OBS"]
    jd = Time(str(mjd), format='mjd', scale='utc').jd
    # hdu1: 1024 1D spectrum
    # hdu2: image 1024 x 512
    # hdu3: image 1024 x 512

    # Structure of the table depends on the object type (STD or not) 
    if obj.strip() == "STD":
        # hdu1
        # 0: wavelength [m]
        # 1: spec model PH        [J*radian/m^3/s]
        # 2: spec model XC        [J*radian/m^3/s] 
        # 3: spec sky             [ADU/s]
        # 4: spec extracted       [ADU/s]
        # 5: spec extracted error [ADU/s]
        # 6: standard star model  [mJy]
        # 7: sensitivity          [mJy]
        data = hdu1.data
        arr = np.array(data, dtype=[
            ("w_m", "f8"), ("specmodel_PH", "f8"), ("specmodel_XC", "f8"),
            ("flux_sky", "f8"), ("flux", "f8"), ("fluxerr", "f8"),
            ("flux_model", "f8"), ("sensitivity", "f8"), 
            ])
    else:
        # hdu1
        # 0: wavelength [m]
        # 1: spec model PH        [J*radian/m^3/s]
        # 2: spec model XC        [J*radian/m^3/s] 
        # 3: spec sky             [ADU/s]
        # 4: spec extracted       [ADU/s]
        # 5: spec extracted error [ADU/s]
        # 6: calibrated flux      [mJy] (not Jy, VISIR manual 1.11 is wrong)
        # 7: flux err             [mJy] (not Jy, VISIR manual 1.11 is wrong)
        data = hdu1.data
        arr = np.array(data, dtype=[
            ("wavelength", "f8"), ("specmodel_PH", "f8"), ("specmodel_XC", "f8"),
            ("flux_sky", "f8"), ("flux_adu", "f8"), ("fluxerr_adu", "f8"),
            ("flux", "f8"), ("fluxerr", "f8"), 
            ])

    df = pd.DataFrame(arr)

    # Remove NaN
    N0 = len(df)
    mask = ~np.isnan(df["flux"])
    df = df[mask]
    N1 = len(df)
    if N0 != N1:
        print(f"  N={N0-N1} data with NaN removed.")

    ## Convert [m] to [micron]
    df["wavelength"] = df["wavelength"]*1e6

    # TODO:Consider light-time correction, exposure time, etc.
    df["jd"] = jd
    df["cflag"] = 0
    df["memo"] = f"VLT_VISIR_{obj}"
    df["code"] = 309
    return df


def plot_spectrum_full(
        df_obj, df_S, df_phot=None, w_fit_range=(7, 13), fit_degree=3, out=None):
    """Plot extracted spectrum with additional information.

    Parameters
    ----------
    df_obj : pandas.DataFrame
        dataframe of object of interest
    df_S : pandas.DataFrame
        dataframe of standard star
    df_phot : pandas.DataFrame, optional
        photometric results of object of interest
    w_fit_range : array-like
        wavelength to be used for the fitting
    fit_degree : int
        degree of polynomial
    out : str, optional
        output filename
    """
    columns = [
        ("flux", "Calibrated Flux", "black", "-"),   
        ("flux_sky", "Sky Background", "blue", "-"),    
        ("dummy", "", "gray", "-"), 
    ]

    def add_telluric_shading(ax):
        for w_min, w_max, label in telluric_features:
            ax.axvspan(w_min, w_max, color='gray', alpha=0.3, label=label)

    # 0: flux
    # 1: sky
    # 2: corrected flux & emissivity
    nrows = 3  
    fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(14, 7), sharex=True)

    dfs = [df_obj, df_S]
    ylabels = ["ADU/s", "Sky", "Corrected Flux"]

    for col_idx, df in enumerate(dfs):
        for row_idx, (col, title, color, ls) in enumerate(columns):
            if row_idx > 2:
                continue  # row_idx: 0=flux, 1=sky, 2=corrected
            ax = axes[row_idx, col_idx]
            if isinstance(col, tuple):
                for c, t, clr, style in zip(col, title, color, ls):
                    if c == "flux_model" or c not in df.columns:
                        continue
                    ax.plot(df["w"], df[c], label=t, lw=1, color=clr, linestyle=style)
            else:
                if col == "flux_model" or col not in df.columns:
                    continue
                ax.plot(df["w"], df[col], label=title, lw=1, color=color, linestyle=ls)
            ax.set_ylabel(ylabels[row_idx], fontsize=9)
            ax.tick_params(axis='both', labelsize=8)
            ax.legend(loc="upper right", fontsize=7)
            ax.grid(True)
            add_telluric_shading(ax)

    # transmissivity & corrected flux
    with np.errstate(divide='ignore', invalid='ignore'):
        transmissivity = np.where(df_S["flux_model"] != 0, df_S["flux"] / df_S["flux_model"], np.nan)
        corrected_flux = np.where(transmissivity != 0, df_obj["flux"] / transmissivity, np.nan)

    # corrected flux plot
    row_corr = 2
    ax_corr = axes[row_corr, 0]
    ax_corr.plot(df_obj["w"], corrected_flux, color="darkgreen", lw=1, zorder=100, label="Corrected Flux (calc)")

    if "flux_cor" in df_obj.columns:
        if "fluxerr_cor" in df_obj.columns:
            ax_corr.errorbar(
                df_obj["w"], df_obj["flux_cor"], yerr=df_obj["fluxerr_cor"],
                fmt="o", markersize=2, color="orange", label="Corrected Flux (from col)",
                alpha=0.7, capsize=1)
        else:
            ax_corr.plot(
                df_obj["w"], df_obj["flux_cor"], color="orange", lw=0.2,
                label="Corrected Flux (from col)")
            
    # Add photometry
    if df_phot is not None:
        label = "Photometry"
        ax_corr.scatter(
            df_phot["wavelength"], df_phot["flux"]*1e3, color="red", label=label)

    ax_corr.set_ylabel("Calibrated Flux [mJy]", fontsize=9)
    ax_corr.set_xlabel("Wavelength [micron]", fontsize=9)
    ax_corr.tick_params(labelsize=8)
    ax_corr.legend(loc="lower right", fontsize=8)
    ax_corr.grid(True)
    add_telluric_shading(ax_corr)

    valid_vals = corrected_flux[~np.isnan(corrected_flux)]
    if len(valid_vals) > 0:
        _, high = np.percentile(valid_vals, [0, 99])
        ax_corr.set_ylim(0, high)

    # fitting
    mask = ~np.isnan(corrected_flux)
    xeff = df_obj["w"][mask]
    yeff = corrected_flux[mask]

    w_min, w_max = w_fit_range
    fit_mask = (xeff >= w_min) & (xeff <= w_max)
    xfit_eff = xeff[fit_mask]
    yfit_eff = yeff[fit_mask]

    polymodel = Polynomial1D(degree=fit_degree)
    linfitter = LinearLSQFitter()
    y_fit_model = linfitter(polymodel, xfit_eff, yfit_eff)

    x_fit = np.linspace(w_min, w_max, 500)
    y_fit = y_fit_model(x_fit)
    ax_corr.plot(x_fit, y_fit, color="red", lw=2, label=f"Poly fit deg={fit_degree}")

    # emissivity
    with np.errstate(divide='ignore', invalid='ignore'):
        emissivity = np.where(y_fit_model(xeff) != 0, yeff / y_fit_model(xeff), np.nan)

    ax_emis = axes[row_corr, 1]
    ax_emis.plot(xeff, emissivity, color="blue", lw=1, label="Emissivity")
    ax_emis.set_ylabel("Emissivity", fontsize=9)
    ax_emis.set_xlabel("Wavelength [micron]", fontsize=9)
    ax_emis.tick_params(labelsize=8)
    ax_emis.legend(loc="upper right", fontsize=8)
    ax_emis.grid(True)
    ax_emis.set_ylim(0.9, 1.1)
    add_telluric_shading(ax_emis)

    plt.tight_layout()
    plt.show(block=False)
    ans = input("Save figure? (y/n): ").strip().lower()
    if ans != 'y':
        plt.close()
    try:
        plt.savefig(out)

    except ValueError:
        print("Not saved. Exiting.")
    plt.close()


def plot_spectrum(
        df, df_phot=None, w_fit_range=(7, 13), fit_degree=3, 
        out_fig=None, out_res=None):
    """Plot extracted spectrum with additional information.

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe of object of interest
    df_phot : pandas.DataFrame, optional
        photometric results of object of interest
    w_fit_range : array-like
        wavelength to be used for the fitting
    fit_degree : int
        degree of polynomial
    out : str, optional
        output filename
    """
    def add_telluric_shading(ax):
        for w_min, w_max, label in telluric_features:
            ax.axvspan(w_min, w_max, color='gray', alpha=0.3, label=label)

    nrows = 1
    # 0: corrected flux & emissivity
    fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(12, 6), sharex=True)

    ax_flux = axes[0]
    ax_emis = axes[1]

    ax_flux.set_ylabel("Flux density [Jy]", fontsize=9)
    ax_flux.set_xlabel("Wavelength [micron]", fontsize=9)
    ax_flux.tick_params(labelsize=8)
    ax_flux.legend(loc="lower right", fontsize=8)
    ax_flux.grid(True)
    add_telluric_shading(ax_flux)

    # From mJy to Jy
    df["flux"] /= 1e3
    df["fluxerr"] /= 1e3

    # Plot observed spectrum        
    ## Add photometry
    if df_phot is not None:
        memo = df_phot["memo"][0]
        label = "Photometry"
        df_phot["width"] = df_phot["wavelength_max"] - df_phot["wavelength_min"]
        ax_flux.errorbar(
            df_phot["wavelength"], df_phot["flux"], 
            xerr=df_phot["width"],
            yerr=df_phot["fluxerr"],
            fmt="o",
            color="red", zorder=100, label=label)

        # Shift spectrum to match photometry
        # TODO: Update
        # Linear interpolation

        w0_phot = np.min(df_phot["wavelength"])
        w1_phot = np.max(df_phot["wavelength"])
        f0_phot = np.mean(df_phot[df_phot["wavelength"]==w0_phot]["flux"])
        f1_phot = np.mean(df_phot[df_phot["wavelength"]==w1_phot]["flux"])
        
        # Flux density of spectrum at w0 and w1
        idx_w0 = (df["wavelength"] - w0_phot).abs().idxmin()
        idx_w1 = (df["wavelength"] - w1_phot).abs().idxmin()
        f0_spec = df.loc[idx_w0, "flux"]
        f1_spec = df.loc[idx_w1, "flux"]

        ratio0 = f0_phot/f0_spec
        ratio1 = f1_phot/f1_spec
        ratio_shift = (ratio0 + ratio1)/2.
        print(f"ratio0, ratio1 = {ratio0:.2f}, {ratio1:.2f}")
        print(f"-> Use average ratio, {ratio_shift:.2f}, to shift")

        # Use first measurement to shift spectrum
        df["flux"] = [y*ratio_shift for y in df["flux"]]
        df["fluxerr"] = [y*ratio_shift for y in df["fluxerr"]]
        
        # Shift using linear interpolation ====================================
        ## If ratio0 == ratio1, it is easy. 
        ## Make a linear model for correction.
        #slope = (ratio1 - ratio0) / (w1_phot - w0_phot)
        #intercept = ratio0 - slope * w0_phot
        #def ratio(w):
        #    return slope * w + intercept

        ## Do correction
        ### Obs
        #ratios = ratio(df["w"])
        #df["flux_cor"] = [y*r for (y, r) in zip(df["flux_cor"], ratios)]
        #df["fluxerr_cor"] = [y*r for (y, r) in zip(df["fluxerr_cor"], ratios)]
        # Shift using linear interpolation ====================================

        label_spec = (
            f"Observed spectrum {memo}\n  "
            f"(shifted to fit photometry by multiplying by {ratio_shift:.2f})")
    else:
        label_spec = "Observed spectrum (not shifted)"
    ## Plot
    ax_flux.errorbar(
        df["wavelength"], df["flux"], yerr=df["fluxerr"],
        fmt="o", markersize=2, color="black", label=label_spec,
        alpha=0.7, capsize=1)

    # Fitting curve
    ## Make fitting curve
    xeff = df["wavelength"]
    yeff = df["flux"]

    w_min, w_max = w_fit_range
    fit_mask = (xeff >= w_min) & (xeff <= w_max)
    xfit_eff = xeff[fit_mask]
    yfit_eff = yeff[fit_mask]

    polymodel = Polynomial1D(degree=fit_degree)
    linfitter = LinearLSQFitter()
    y_fit_model = linfitter(polymodel, xfit_eff, yfit_eff)
     
    # Use the same wavelengths as observations
    y_fit = y_fit_model(xeff)

    # Plot
    label_fit = f"Poly fit deg={fit_degree}"
    ax_flux.plot(
        xeff, y_fit, color="gray", ls="dashed", lw=2, label=label_fit)
    
    # Change plot range
    #valid_vals = corrected_flux[~np.isnan(corrected_flux)]
    #if len(valid_vals) > 0:
    _, high = np.percentile(df["flux"], [0, 99])
    #ax_flux.set_ylim(0, high)
    ax_flux.set_ylim(0, 100)
    ax_flux.set_xlim(7, 14)
    
    # emissivity
    with np.errstate(divide='ignore', invalid='ignore'):
        emissivity = np.where(y_fit_model(xeff) != 0, yeff / y_fit_model(xeff), np.nan)
    df["emissivity"] = emissivity
    df[f"flux_poly{fit_degree}"] = y_fit
    df["scalefactor"] = ratio_shift

    add_telluric_shading(ax_emis)
    ax_emis.plot(xeff, emissivity, color="black", lw=1, label="Emissivity")
    ax_emis.set_ylabel("Emissivity", fontsize=9)
    ax_emis.set_xlabel("Wavelength [micron]", fontsize=9)
    ax_emis.tick_params(labelsize=8)
    ax_emis.grid(True)
    ax_emis.set_ylim(0.90, 1.10)
    ax_emis.legend(fontsize=12).get_frame().set_alpha(1.0)
    ax_flux.legend(fontsize=12).get_frame().set_alpha(1.0)

    plt.tight_layout()
    plt.show(block=False)
    ans = input("Save figure? (y/n): ").strip().lower()
    if ans != 'y':
        plt.close()
    try:
        plt.savefig(out_fig)
        print(f"  Figure is saved as {out_fig}")
        # Save spectroscopy
        col_save = ["jd", "wavelength", "flux", "fluxerr", "scalefactor", "code", "cflag", "memo"]
        
        df = df[col_save]
        df.to_csv(out_res, sep=" ", index=False)
        print(f"  Result is saved as {out_res}")

    except ValueError:
        print("Not saved. Exiting.")
    plt.close()


if __name__ == "__main__":
    parser = ap(description="Extract spectrum.")
    parser.add_argument(
        "fi", type=str,
        help="Fits file")
    parser.add_argument(
        "--fi_S", type=str, default=None,
        help="Fits file of standard star (necessary with full option)")
    parser.add_argument(
        "--f_phot", type=str, default=None,
        help="Result of photometry")
    parser.add_argument(
        "--full", action="store_true", default=False,
        help="Plot full results.")
    parser.add_argument(
        "--out_fig", type=str, default="spectrum.pdf",
        help="Output figure name")
    parser.add_argument(
        "--out_res", type=str, default="spectrum.txt",
        help="Output filename")
    args = parser.parse_args()
    
    # Extract
    df = extract_1dspec(args.fi)
    

    # Photometry if exists
    if args.f_phot:
        df_phot = pd.read_csv(args.f_phot, sep=" ")
    else:
        df_phot = None
    
    if args.full:
        df_S = extract_1dspec(args.fi_S)
        plot_spectrum_full(
            df, df_S, df_phot=df_phot, w_fit_range=(8.2, 13), fit_degree=3, out=args.out)
    else:
        plot_spectrum(
            df, df_phot=df_phot, w_fit_range=(8.2, 13), 
            fit_degree=3, out_fig=args.out_fig, out_res=args.out_res)
