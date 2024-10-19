import pandas as pd


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
