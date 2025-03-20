import matplotlib.pyplot as plt



mycolor = [
    # red
    "#AD002D", 
    # blue
    "#1e50a2", 
    # green
    "#69821b", 
    # orange
    #"#ec6800", 
    # purple
    "#f055f0", 
    # gray
    "#afafb0", 
    "#0095b9", 
    "#89c3eb"] 

# linestyle
myls = ["solid", "dashed", "dashdot", "dotted", (0, (5, 3, 1, 3, 1, 3)), 
        (0, (4,2,1,2,1,2,1,2))]
# marker
mymark = ["o", "^", "s", "D", "*", "v", "<", ">", "h"]

def plot_growth(df, t_lim, snr_lim, out):
    """
    Plot growth curve

    Parameters
    ----------
    df : pandas.DataFrame
        photometric results of standard star
    t_lim : float
        Effective exposure time in s to calculate limiting flux
    snr_lim : float
        SNR of limiting flux
    """

    fig = plt.figure(figsize=(12, 8))
    ax_u = fig.add_axes([0.15, 0.6, 0.7, 0.35])
    
    # Flux in aperture
    ax_u.set_ylabel("Flux [ADU]")
    ax_u.scatter(df["radius"], df["flux"], color="black", label="Flux in aperture")
    
    # SNR with red
    ax_SNR = ax_u.twinx()
    ax_SNR.set_ylabel("SNR", color="red")
    ax_SNR.spines['right'].set_color('red')
    ax_SNR.tick_params(axis='y', colors='red')
    ax_SNR.scatter(df["radius"], df["snr"], color="red", marker="+", label="SNR")
  
    
    ax_l = fig.add_axes([0.15, 0.15, 0.7, 0.35])
    ax_l.set_xlabel("Aperture radius [pix]")
    band = df["band"][0]
    ax_l.set_ylabel(f"{band} limiting flux [mJy]")
    label_lim = f"Flux to achieve SNR of {snr_lim} with exposure of {t_lim} s"
    ax_l.scatter(df["radius"], df["fluxlim"], color="black", label=label_lim)
    ax_SNR.legend(loc="upper left") 
    ax_u.legend(loc="lower right")
    ax_l.legend(loc="upper right")
    
    x, y = -0.12, 0.5
    ax_u.yaxis.set_label_coords(x, y)
    ax_l.yaxis.set_label_coords(x, y)

    plt.savefig(out)
    plt.close()
