#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Show image and save it.
"""
from argparse import ArgumentParser as ap
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt


if __name__ == "__main__":
    parser = ap(description="Make a stacked image for VLT/VISIR.")
    parser.add_argument(
        "fi", type=str,
        help="Fits to be analyzed.")
    parser.add_argument(
        "--cmap", default="inferno",
        help="Color map")
    parser.add_argument(
        "--sigma", type=int, default=3,
        help="Dynamic range of the image")
    parser.add_argument(
        "--out_image", type=str, default="image_visir.jpg",
        help="output image")
    args = parser.parse_args()
    
    # # Read 3 colors images ====================================================
    # img_list = []
    # for fi in args.fi:
    #     hdu = fits.open(fi)
    #     img = hdu[0].data
    #     hdr = hdu[0].header
    #     img = img.byteswap().newbyteorder()
    #     # Normalize
    #     img_list.append(img)
    # rows, cols = img_list[0].shape
    # # 3 colors with same dimension
    # result = np.zeros((rows, cols, 3))
    # result[..., 0] = img_list[0]
    # result[..., 1] = img_list[1]
    # result[..., 2] = (img_list[0] + img_list[1])/2.
    # ax.imshow(result)
    # # Read 3 colors images ====================================================
    
    hdu = fits.open(args.fi)
    img = hdu[0].data
    hdr = hdu[0].header
    img = img.byteswap().newbyteorder()

    ## Dynamic range
    sigma = args.sigma
    stdsub = np.std(img)
    vmin, vmax = -sigma*stdsub, sigma*stdsub
    
    # Plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    ax.set_xlabel("x [pix]")
    ax.set_ylabel("y [pix]")
    ax.imshow(img, cmap=args.cmap, vmin=vmin, vmax=vmax)

    plt.savefig(args.out_image)
