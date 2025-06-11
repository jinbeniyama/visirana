#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Plot spectrum.
"""
import os 
import sys
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  

from myplot import mycolor, myls

if __name__ == "__main__":
    parser = ap(description="Plot spectrum.")
    parser.add_argument(
        "sp", nargs="*", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "--label", nargs="*", type=str, default=None,
        help="Label of spectrum")
    parser.add_argument(
        "--color", nargs="*", type=str, default=None,
        help="Color of spectrum")
    parser.add_argument(
        "--ls", nargs="*", type=int, default=None,
        help="'Index' of line style of spectrum")
    parser.add_argument(
        "--out", type=str, default="spectra.jpg",
        help="Output filename")
    args = parser.parse_args()
    

    fig = plt.figure(figsize=(12, 4))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    
    for idx_sp, sp in enumerate(args.sp):
        df = pd.read_csv(sp, sep=" ")

        if args.label is None:
            label = None
        else:
            label = args.label[idx_sp]

        if args.color is None:
            col = mycolor[idx_sp]
        else:
            col = args.color[idx_sp]

        if args.ls is None:
            ls = myls[idx_sp]
        else:
            ls = myls[args.ls[idx_sp]]

        ax.plot(df["x"], df["flux"], label=label, color=col, ls=ls)
    
    ax.set_ylabel("Flux")
    ax.set_xlabel("x [pixel]")
    ax.legend()
    plt.show(block=False)

    ans = input("Save figure? (y/n): ").strip().lower()
    if ans != 'y':
        plt.close()
    try:
        plt.savefig(args.out)

    except ValueError:
        print("Not saved. Exiting.")
    plt.close()
