#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Stack raw images taken with VLT/VISIR.

ABAB or ABBA nodding can be specified.
"""
from argparse import ArgumentParser as ap
import os
import numpy as np
import astropy.io.fits as fits


def main(args):
    """This is the main function called by the `nodstack.py` script.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments passed from the command-line as defined below.
    """

    sigma = 3
    Vshift = None
    if args.f_list:
        fA_list, fB_list = [], []
        for idx, f in enumerate(args.f_list):
            # A-B-B-A nodding 
            if args.nodtype == "ABBA":
                if ((idx % 4) == 0) or (idx % 4 == 3): 
                    fA_list.append(f)
                else:
                    fB_list.append(f)
            # A-B-A-B nodding 
            if args.nodtype == "ABAB":
                if ((idx % 2) == 0): 
                    fA_list.append(f)
                else:
                    fB_list.append(f)
    else:
        fA_list, fB_list = args.fA_list, args.fB_list

    assert len(fA_list)==len(fB_list), "Check inputs."

    fAstack_list, fBstack_list, substack_list = [], [], []
    for idx_cycle, (fA, fB) in enumerate(zip(fA_list, fB_list)):
        
        hdu_A = fits.open(fA)
        hdu_B = fits.open(fB)
        N_hdu = len(hdu_A)
        # For old data, N_hdu == 1.
        # The hdu has 3-d image.
        # For new data, N_hdu == 4.
        # The hdu0 has header wo/ image.
        # The hdu1 has a 2-d image (chop1).
        # The hdu2 has a 2-d image (chop2).
        # The hdu3 has a 2-d image (chop1 - chop2).

        if len(hdu_A) == 1:
            if idx_cycle == 0:
                print("Analyze old data w/stacking.")
            hdr_A = hdu_A[0].header
            image_A = hdu_A[0].data

            hdr_B = hdu_B[0].header
            image_B = hdu_B[0].data

            # Stack chopping pair 1 and 2
            image_A_c1 = np.sum(image_A[0::2,::,::], axis=0)
            image_A_c2 = np.sum(image_A[1::2,::,::], axis=0)
            # Subtract 
            image_A_c12 = image_A_c1 = image_A_c2

            # Stack chopping pair 1 and 2
            image_B_c1 = np.sum(image_B[0::2,::,::], axis=0)
            image_B_c2 = np.sum(image_B[1::2,::,::], axis=0)
            # Subtract 
            image_B_c12 = image_B_c1 = image_B_c2

            ## Shift
            #if args.vshift:
            #    vx, vy = vshift
            #    image_A = np.roll(image_A, vx*idx_cycle, axis=0)
            #    image_A = np.roll(image_A, vy*idx_cycle, axis=1)
            #    image_B = np.roll(image_B, vx*idx_cycle, axis=0)
            #    image_B = np.roll(image_B, vy*idx_cycle, axis=1)

        else:
            if idx_cycle == 0:
                print("Analyze new data wo/stacking.")
            hdr_A = hdu_A[0].header
            image_A_c12 = hdu_A[3].data
            if args.cut:
                image_A_c12 = cut4VISIR(image_A_c12)


            hdr_B = hdu_B[0].header
            image_B_c12 = hdu_B[3].data
            if args.cut:
                image_B_c12 = cut4VISIR(image_B_c12)

        # Subtract
        image_sub = image_A_c12 - image_B_c12

        # Save 
        #fAstack_list.append(image_A_c12)
        #fBstack_list.append(image_B_c12)
        substack_list.append(image_sub)

    # Stacking 
    #imageA = np.sum(fAstack_list, axis=0)
    #imageB = np.sum(fAstack_list, axis=0)    
    imagesub = np.sum(substack_list, axis=0)

    # Save
    fits0 = fA_list[0]
    hdu0 = fits.open(fits0)
    hdr0 = hdu0[0].header
    sta = fits.PrimaryHDU(data=imagesub, header=hdr0)

    # Add history
    hdr = sta.header
    hdr.add_history(
        f"[nodstack] header info. is inherited from {fits0}")
    for idx, f in enumerate(fA_list):
        hdr.add_history(
            f"[nodstack] fits {idx+1} (nodding A) : {f}")
    for idx, f in enumerate(fB_list):
        hdr.add_history(
            f"[nodstack] fits {idx+1} (nodding B) : {f}")
    # Card 'ESO DET CHIP PXSPACE' is not FITS standard
    key_del = "ESO DET CHIP PXSPACE"
    if key_del in hdr:
        del hdr[key_del]
    
    # Check the keyword for GAIN
    kwd_gain = "HIERARCH ESO DET CHIP GAIN"
    if len(hdu_A) == 1:
        if kwd_gain in hdr:
            pass
        else:
            print("WARNING: No GAIN keyword was found.")
            print("WARNING: Set GAIN as 1.")
            hdr[kwd_gain] = 1.0
    else:
        # For new fits
        # Add gain from hdr1
        hdr1 = hdu0[1].header
        gain = hdr1[kwd_gain]
        hdr[kwd_gain] = gain
    sta.writeto(args.out, overwrite=True)


if __name__ == "__main__":
    parser = ap(description="Make a stacked image for VLT/VISIR.")
    parser.add_argument(
        "nodtype", type=str,
        help="ABBA or ABAB")
    parser.add_argument(
        "--f_list", type=str, nargs="*", default=None,
        help="Images taken in ABBA or ABAB sequence")
    parser.add_argument(
        "--fA_list", type=str, nargs="*", default=None,
        help="Images taken on nodding position A")
    parser.add_argument(
        "--fB_list", type=str, nargs="*", default=None,
        help="Images taken on nodding position B")
    parser.add_argument(
        "--cut", action="store_true", default=False,
        help="Cut useless pixels")
    parser.add_argument(
        "--vshift", type=float, default=None,
        help="Velocity of shift in pix/frame")
    parser.add_argument(
        "--out", type=str, default="stacksub.fits",
        help="output fits file")
    args = parser.parse_args()
    
    main(args)
