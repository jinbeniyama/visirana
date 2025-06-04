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
from scipy.interpolate import interp1d


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


def extract_sky_template(fi):
    """Extract atmospheric emission spectrum.

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
    # 1: emission
    data = hdu1.data

    arr = np.array(data, dtype=[("w_m", "f8"), ("emi", "f8")])
    df = pd.DataFrame(arr)
    return df


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
    parser = ap(description="Plot model flux of standard star.")
    parser.add_argument(
        "fi", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "fi_sky", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "fi_qe", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "obj", nargs="*", type=str,
        help="Objects of interest.")
    parser.add_argument(
        "--xr", nargs=2, type=float, default=None,
        help="X range")
    parser.add_argument(
        "--out", type=str, default="modelflux.jpg",
        help="Output filename")
    args = parser.parse_args()
    
    # Extract
    df_sky = extract_sky_template(args.fi_sky)
    ## Convert [m] to [micron]
    df_sky["w"] = df_sky["w_m"]*1e6
    interp_sky = interp1d(
        df_sky['w'], df_sky['emi'], kind='linear', 
        bounds_error=False, fill_value=np.nan)

    # Extract
    df_qe = extract_qe_template(args.fi_qe)
    ## Convert [m] to [micron]
    df_qe["w"] = df_qe["w_m"]*1e6
    interp_qe = interp1d(
        df_qe['w'], df_qe['eff'], kind='linear', 
        bounds_error=False, fill_value=np.nan)

    # ???
    #df["emi"] = interp_sky(df["w"])
    #df["eff"] = interp_qe(df["w"])
    #df["flux_jy_cor"] = df["flux_jy"] * df["emi"] * df["eff"]
    #ax.plot(df["w"], df["flux_jy_cor"], label=f"Model x transmission {ob}")

    # 3行1列のプロット領域を用意
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 10), sharex=True)
    
    # 軸の参照
    ax_flux, ax_qe, ax_sky = axes
    
    # 上段: 各 obj のモデルフラックスを Jy 単位でプロット
    for ob in args.obj:
        df = extract_modelflux(args.fi, ob)
        if df is None:
            print(f"{ob} not found in {args.fi}.")
            sys.exit()
    
        # 単位変換
        df["w_m"] = df["w"] * 1e-6  # [micron] → [m]
        c = 2.998e8  # [m/s]
        df["flux_jy"] = (df["w_m"]**2 / c * df["flux"]) / 1e-26  # [Jy]
    
        ax_flux.plot(df["w"], df["flux_jy"], label=f"Model {ob}")
    
    ax_flux.set_ylabel("Flux density [Jy]", fontsize=9)
    ax_flux.legend(fontsize=8)
    ax_flux.grid(True)
    ax_flux.tick_params(labelsize=8)
    
    # 中段: df_qe のプロット
    ax_qe.plot(df_qe["w"], df_qe["eff"], color='blue', label="QE")
    ax_qe.set_ylabel("Efficiency", fontsize=9)
    ax_qe.legend(fontsize=8)
    ax_qe.grid(True)
    ax_qe.tick_params(labelsize=8)
    
    # 下段: df_sky のプロット
    ax_sky.plot(df_sky["w"], df_sky["emi"], color='gray', label="Sky Emission")
    ax_sky.set_ylabel("Emission", fontsize=9)
    ax_sky.set_xlabel("Wavelength [micron]", fontsize=9)
    ax_sky.legend(fontsize=8)
    ax_sky.grid(True)
    ax_sky.tick_params(labelsize=8)

    if args.xr is None:
        pass
    else:
        x0, x1 = args.xr
        for ax in fig.axes:
            ax.set_xlim([x0, x1])
    
    # 全体調整
    plt.tight_layout()
    plt.show(block=False)


    #fig = plt.figure()
    #ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    #ax.set_xlabel("Wavelength [micron]")
    #ax.set_ylabel("Flux density [Jy]")
    #for ob in args.obj:
    #    df = extract_modelflux(args.fi, ob)
    #    if df is None:
    #        print(f"{ob} not found in {args.fi}.")
    #        sys.exit()
    #    # Original unit is [micron] and [W/m^2/m] (See 7.2 Spectroscopy Standard Stars of visir-pipeline-manual-1.11.pdf)
    #    ## Convert [micron] to [m]
    #    df["w_m"] = df["w"]*1e-6
    #    ## Convert [W/m^2/m] to [Jy]
    #    # speed of light [m/s]
    #    c = 2.998e8
    #    df["flux_jy"] = (df["w_m"]**2/c * df["flux"]) / 1e-26
    #    ax.plot(df["w"], df["flux_jy"], label=f"Model {ob}")
    #ax.legend()
    #plt.show(block=False)



    ans = input("Save figure? (y/n): ").strip().lower()
    if ans != 'y':
        plt.close()
    try:
        plt.savefig(args.out)

    except ValueError:
        print("Not saved. Exiting.")
    plt.close()
