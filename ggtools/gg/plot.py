from os import path,makedirs
import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def plot_at_northpole(lons,lats,data,fig_name,ylabel=None,magnify=None):
    '''
    Plot grid data at the northpole. Note: The value of the grid data was magnified 1000 times during the plot process. 

    Usage:
    plot_at_northpole(lons,lats,data,fig_name)
    plot_at_northpole(lons,lats,data,fig_name,1e3)

    Inputs:
    lons -> [float array] logitudes
    lats -> [float array] latitudes
    data -> [float 2d array] grid data
    fig_name -> [str] figure name
    ylabel -> [str] ylabel, such as '$10^{-3}$ [mm w.e.]'

    Parameters:
    magnify -> [optional, float, default = None] If None, 1 is taken.
    
    Outputs: A png(200dpi) image stored in the figures directory 
    '''
    # set a directory to store images
    fig_dir = 'figures/'
    if not path.exists(fig_dir): makedirs(fig_dir) 

    plt.clf()
    fig = plt.figure(dpi=200)

    # set the projection to PlateCarree
    proj = ccrs.NearsidePerspective(0,90,2e6)
    ax = fig.add_subplot(1, 1, 1,projection = proj)

    # set the gridlines to dashed line and the transparency to 0.7
    gl = ax.gridlines(xlocs=np.linspace(0,360,7),ylocs=np.linspace(90,0,10),linestyle='--',alpha=0.7)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # add coastline, rivers, lakes
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.RIVERS) 
    ax.add_feature(cfeature.LAKES) 
    
    if magnify is None: magnify = 1
    XX, YY = np.meshgrid(lons,lats)
    Z = data*magnify
    
    # calculate the maximum of the abs(Z)
    abs_Z_max = np.abs(Z).max()

    # set the contour levels
    Z_levels = np.linspace(-abs_Z_max,abs_Z_max, 61)
    CS = ax.contourf(XX,YY,Z,levels = Z_levels,extend='both',cmap=plt.cm.RdBu_r,zorder=0,transform = ccrs.PlateCarree())
    
    cbar = plt.colorbar(CS,extend='both',format='%.0f',shrink=0.9)
    cbar.ax.set_ylabel(ylabel,fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    return plt.savefig(fig_name,bbox_inches='tight')