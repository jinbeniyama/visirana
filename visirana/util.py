import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

import sep

#   wmin, wmax, label
telluric_features = [
    (9.2, 10.1, "Ozone"),
    (6.0, 8.2, "Water?"),
    (13.0, 14.3, "Water?"),
]

def cut4VISIR(image):
    """
    Cut useless pixels for clarity.
    Note: 
    Some useless pixels could exist along y axis as well.

    Parameter
    ---------
    image : 2-d array-like
        input image

    Return
    ------
    image_cut : 2-d array-like
        output image
    """
    image_cut = image[:, 25:875]
    return image_cut


def Perr(x, gain):
    """
    Calculate Poission error.

    Parameters
    ----------
    x : float
        count in ADU
    gain : float
        inverse gain e/ADU

    Return
    ------
    perr : float
        Poisson erorr in ADU
    """
    # x [ADU] -> x*gain [e]
    # root(x*gain) [e] -> root(x*gain)/gain [ADU]
    perr = (x*gain)**0.5/gain
    return perr


def get_VISIR_standard():
    """
    Read standard stars
    TODO: Add other fluxes rather than B8.7, B10.7, and B11.7.
    """
    f = "/Users/beniyama/research/1998KY26_VLT/data/zerop_cohen_Jy.txt"
    with open(f, "r") as f:
        obj_list, B8_7_list, B10_7_list, B11_7_list = [], [], [], []
        lines = f.readlines()
        for l in lines[1:]:
            obj   = l[0:11]
            B8_7  = l[258:258+9]
            B10_7 = l[278:278+9]
            B11_7 = l[288:288+9]
    
            obj_list.append(obj.replace(" ", ""))
            B8_7_list.append(float(B8_7))
            B10_7_list.append(float(B10_7))
            B11_7_list.append(float(B11_7))
    
    df_stan = pd.DataFrame({
        "obj":obj_list,
        "B8.7":B8_7_list,
        "B10.7":B10_7_list,
        "B11.7":B11_7_list,
        })
    return df_stan


def obtain_winpos(data, x, y, radius, nx, ny):
    """ 
    Obtain windowed centroid xwin and ywin.
    Note: x and y should be numpy.ndarray

    Parameters
    ----------
    data : numpy.ndarray
      2-d image data
    x, y : numpy.ndarray
      location(s) of object(s)
    radius : float
      scale related to the area where searched the centroid
    """

    # Radius should be large enough to obtain total flux
    # wpos_param: constant to convert `0.5*FWHM` to `sigma` (FWHM = 2.35*sigma) 
    # 0.5*FWHM*wpos_param = 0.5*FWHM*2/2.35 = sigma
    wpos_param  = 2.0/2.35
    # Half flux radius
    frad_frac   = 0.5
    frad_subpix = 5
    frad_ratio  = 5.0

    # Do photometry to obtain all flux
    # Note1: err and gain are needless for flux estimation
    # Note2: Sky background should be subtracted (?) 

    # Must need
    flux,fluxerr,eflag = sep.sum_circle(data, x, y, r=radius)

    # Use only objects with nonzero eflag (?)
    # Create array like [radius, radius, ... , radius]
    # Note: radius is not used in flux_radius when normflux is used (?)
    # r : flux radius, i.e., `0.5*fwhm!`
    radius = np.full_like(x, radius)
    r, flag = sep.flux_radius(
        data, x, y, radius, frad_frac, normflux=flux, subpix=frad_subpix)
    # r (0.5*FWHM) to sigma (see above)
    sigma      = wpos_param*r
    sigma_mean = np.mean(sigma)
    sigma_std  = np.std(sigma)
    #print(f"  N={len(sigma)}, calculated in obtain_winpos")
    #print(f"  sigma: {sigma_mean:.1f}+-{sigma_std:.1f}")
    #print(f"  FWHM : {2.35*sigma_mean:.1f} (sigma times 2.35)")

    # wflag is always 0 if mask=None
    # Search winpos with estimated sigma

    # Search narrow region
    # sigma = 0.5*sigma
    # Dramatically works bad for faint objects
    xwin, ywin, wflag = sep.winpos(data, x, y, sigma)

    
    # If the differences of coordinates are larger than ratio_diff*radius,
    # for objects more than ratio_obj*N_obj,
    # print a warning message.
    # original values as xwin and ywin with eflag_win = 1
    ratio_diff = 0.3
    ratio_obj  = 0.5
    diff = np.sqrt((xwin-x)**2 + (ywin-y)**2)
    # 1 for large diff, 0 for small diff
    flag = np.where(diff > ratio_diff*radius, 1, 0)
    ratio_large_diff= np.sum(flag)/flag.size 
    if ratio_large_diff > ratio_obj:
        print(f"      Large winpos correct ratio detected :{ratio_large_diff:.1f}")
        print(f"      This is just a caution. Please check wcs information etc.")

    # Insert original value when xwin and ywin are outside of FoV
    xwin = [x if (x < nx) and (x > 0) and (y < ny) and (y > 0) else x0 for x,y,x0 in zip(xwin,ywin,x)]
    ywin = [y if (x < nx) and (x > 0) and (y < ny) and (y > 0) else y0 for x,y,y0 in zip(xwin,ywin,y)]
    return xwin, ywin, flag


def remove_background2d(image, mask, bw=64, fw=3):
    """
    Remove background from 2D FITS

    Parameters
    ----------
    image : array-like
        input image to be background subtracted
    mask : array-like
        mask array 
    bw : int
        box width 
    fw : int
        filter width 

    Returns
    -------
    image : array like
        background subtracted image
    bg_info : dict
        background info.
    """
    bg_engine = sep.Background(
        image, mask=mask, bw=bw, bh=bw, fw=fw, fh=fw)
    bg_engine.subfrom(image)
    bg_global = bg_engine.globalback
    bg_rms = bg_engine.globalrms
    bg_info = {'level': bg_global, 'rms': bg_rms}
    return image, bg_info


def moffat_fit_oneside(x, a, sigma, beta):
    """
    Moffat fit function for scipy.optimize. (oneside)
    y = a*(1+x**2/sigma**2)**(-beta) (x >= 0)
    """
    return a*(1+x**2/sigma**2)**(-beta)*(x>=0) + 0*(x<0)


def gauss_fit_oneside(x, a, sigma):
    """
    Gaussian fit function for scipy.optimize. (oneside)
    y = a*exp(-0.5*x**2/sigma**2)

    """
    return a*np.exp(-0.5*(x-mu)**2/sigma**2)


def fit_psf(x, y, yerr=None, fittype="Gauss", oneside=False):
    """
    Fit and return parameter of psf.

    Parameters
    ----------
    x, y : array-like
        1-d array for coordinate and pixel value
    yerr : array-like
        1-d array for uncertainty
    fittype : str
        Gauss or Moffat
    oneside : bool
        Set True for radial profile fitting

    Return
    ------
    param : array-like
        fitting parameters
    """
    
    # p0 : initial values

    if fittype == "Gauss":
        if oneside:
            func_fit = gauss_fit_oneside
            p0 = [np.max(y), 1]
        else:
            func_fit = gauss_fit
            p0 = [np.max(y), np.mean(x), 1]
    elif fittype == "Moffat":
        if oneside:
            func_fit = moffat_fit_oneside
            p0 = [np.max(y), 1, 2]
        else:
            func_fit = moffat_fit
            p0 = [np.max(y), np.mean(x), 1, 2]
    
    param, cov = curve_fit(func_fit, x, y, sigma=yerr, p0=p0)
    return param


def radial_profile(image, center, bgerr, gain):
    """
    Create radial profile from the center with certain width.

    Parameters
    ----------
    image : 2d array-like
        2d array for coordinates
    center : tuple
        x and y coordinates
    bgerr : float
        background error in AUD
    gain : float
        inverse gain e/AUD

    Return
    ------
    y : array-like
        radial profile
    yerr : array-like
        error of y
    """

    # Create 2-d grid
    x, y = np.indices((image.shape))
    nx, ny   = image.shape
    c_x, c_y =  center
    print(f"  Shape of cut image: {nx} x {ny}")
    print(f"  Object center ({c_x}, {c_y})")
    # Calculate distance from the center
    r = np.sqrt((x-c_x)**2 + (y-c_y)**2)
    # Convert to integer
    r = r.astype(int)

    # tbin: Total number
    # nr  : number of pixels in the range
    # ex. nr=3 (number of pixel in the range), tbin=4000 (total in the bin),
    #     -> y = 4000/3

    # Input should be integer!
    # Radial profile weighting by pixel values in each 1 pix radius
    #   (r: 1-d radial profile in integer)
    tbin = np.bincount(r.ravel(), weights=image.ravel())
    # nr: number of pixels in each pixel bin (norm factor)
    nr = np.bincount(r.ravel())
    # Normalized by pixels
    y = tbin / nr
    y = y.tolist()
    # Poisson error for electron, not ADU
    # F(ADU) * gain(e/AUD) = F (e)
    # Ferr(e) = (F(e))**0.5
    # Ferr(ADU) = Ferr(e) / gain(e/ADU)
    perr = [np.sqrt(abs(v)*gain)/gain for v in y]
    yerr = [np.sqrt(bgerr**2 + pe**2) for pe in perr]
    return y, yerr


