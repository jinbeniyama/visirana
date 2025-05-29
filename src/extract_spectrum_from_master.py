#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Extract spectrum from M....fits.
"""
import os 
import sys
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  
import astropy.io.fits as fits


def extract_1dspec(fi):
    """Extract calibrated spectrum.

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
    hdu1 = hdu[1]
    # hdu1: 1024 1D spectrum
    # hdu2: image 1024 x 512
    # hdu3: image 1024 x 512
    
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
    df = pd.DataFrame(arr)
    return df


if __name__ == "__main__":
    parser = ap(description="Extract spectrum.")
    parser.add_argument(
        "fi", type=str,
        help="Fits file")
    parser.add_argument(
        "--out", type=str, default="spectrum.txt",
        help="Output filename")
    args = parser.parse_args()
    
    # Extract
    df = extract_1dspec(args.fi)
    
    ## Convert [m] to [micron]
    df["w"] = df["w_m"]*1e6
    

    # プロット対象と設定を定義
    columns = [
        ("flux", "Calibrated Flux", "black", "solid"),
        ("flux_sky", "Sky", "gray", "dashed"),
        ("flux_model", "Flux Model", "blue", "dotted"),
        # 最後のsubplotにまとめて2つ描画
        (("specmodel_PH", "specmodel_XC"), ("Specmodel PH", "Specmodel XC"), ("green", "red"), ("dashdot", "solid")),
    ]
    
    # FigureとAxesの用意（4行1列）
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10, 10), sharex=True)
    axes[0].set_ylabel("ADU/s", fontsize=10)
    axes[1].set_ylabel("ADU/s", fontsize=10)
    axes[2].set_ylabel("mJy", fontsize=10)
    axes[3].set_ylabel("J*radian/m^3/s", fontsize=10)
    
    # 各subplotにプロット
    for i, (col, title, color, ls) in enumerate(columns):
        ax = axes[i]
        if isinstance(col, tuple):
            # specmodel_PH と specmodel_XC を同じsubplotに描画
            for c, t, clr, style in zip(col, title, color, ls):
                ax.plot(df["w"], df[c], label=t, lw=1, color=clr, linestyle=style)
        else:
            ax.plot(df["w"], df[col], label=title, lw=1, color=color, linestyle=ls)
        
        ax.legend(loc="upper right")
        ax.tick_params(axis='both', labelsize=8)
        ax.grid(True)
    
    # x軸ラベルは最後にだけ付ける
    axes[-1].set_xlabel("Wavelength [micron]")
    
    plt.tight_layout()
    plt.show(block=False)


    ans = input("Save figure? (y/n): ").strip().lower()
    if ans != 'y':
        plt.close()
    try:
        plt.savefig(args.out)

    except ValueError:
        print("Not saved. Exiting.")
    plt.close()
