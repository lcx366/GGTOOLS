from os import path,makedirs
import numpy as np
from pyshtools.shclasses import SHGrid
from pyshtools.expand import MakeGridDH,MakeGridPoint

from .filter import filter_ddk,filter_gaussian
from .plot import plot_at_northpole 
import matplotlib.pyplot as plt  

def ddk_gaussian(filter_type,D,visible = None):
    '''
    Given a specific type of DDK filter and the maximum SHC degree number, evaluate the 'equivalent' Gaussian filter radius. 
    Different Gaussian filter radius corresponds to a different correlation between DDK filer and Gaussian filter.
    According to the rule that the largest correlation corresponds to the 'equivalent' radius, this program try to find out and visualize them.

    Usage:
    ddk_gausian('DDK5',96)
    ddk_gausian('DDK3',60,'visible')

    Inputs:
    filter_type -> [str] Types of DDK filter. Available options are 'DDK1' to 'DDK8'.
    D -> [int] Degree and order of SHC

    Parameters:
    visible -> [optional, str, default = None] If None, the visualization of the Point Spreading Function(PSF), DDK filtered PSF, Gaussian filtered PSF at the northpole as well as their cross section will be closed. If 'visible', they will be visualized by outputing an image.
    
    Outputs:
    Print the 'equivalent' Gaussian filter radius and the correlation coefficient between the given DDK filer and the Gaussian filter with the 'equivalent' radius. 
    Also, images of the Point Spreading Function(PSF), DDK filtered PSF, Gaussian filtered PSF at the northpole as well as their cross section can be generated in the figures directory. 
    ''' 
    # create a cap with an angular radius of 0.1 degree at the North Pole
    cap_at_equator_grids_class = SHGrid.from_cap(0.1,0,0,D)
    cap_at_equator_coeffs_class = cap_at_equator_grids_class.expand()
    cap_at_northpole_coeffs_class = cap_at_equator_coeffs_class.rotate(0, 90, 0)
    cap_at_northpole_grids_class = cap_at_northpole_coeffs_class.expand()
    cap_at_northpole_coeffs = cap_at_northpole_coeffs_class.coeffs
    cap_at_northpole_grids = cap_at_northpole_grids_class.data
    
    # calculate the correlations between the DDK filter and Gaussian filter
    filt_SHC_ddk = filter_ddk(filter_type,cap_at_northpole_coeffs)
    ddk_power = np.sum(filt_SHC_ddk**2)
    corrs = [] 

    # set the range of the 'equivalent' filter radius
    rs = range(20,550,5)

    for r in rs:
        filt_SHC_gau = filter_gaussian(r,cap_at_northpole_coeffs)
        gau_power = np.sum(filt_SHC_gau**2)
        ddk_gau_crosspower = np.sum(filt_SHC_ddk*filt_SHC_gau)
        alpha = ddk_gau_crosspower/np.sqrt(ddk_power*gau_power)
        corrs.append(alpha) 
    corrs = np.array(corrs) 
    max_corrs_index = np.argmax(corrs)
    filt_SHC_gau = filter_gaussian(rs[max_corrs_index],cap_at_northpole_coeffs)
    print('Correlation: {:.4f}\nApproximate equivalent Gaussian filter radius for {:s}: {:d}'.format(corrs[max_corrs_index],filter_type,rs[max_corrs_index]))

    thetas = np.arange(0,20,0.1)
    value,value_ddk,value_gau = [],[],[]
    for theta in thetas:
        value.append(MakeGridPoint(cap_at_northpole_coeffs,90-theta,0))
        value_ddk.append(MakeGridPoint(filt_SHC_ddk,90-theta,0))
        value_gau.append(MakeGridPoint(filt_SHC_gau,90-theta,0))
    value = np.array(value)    
    value_ddk = np.array(value_ddk)    
    value_gau = np.array(value_gau)
    
    # plot the Point Spreading Function(PSF), DDK filtered PSF, Gaussian filtered PSF at the northpole as well as their cross section
    if visible is not None:

        fig_dir = 'figures/'
        if not path.exists(fig_dir): makedirs(fig_dir) 

        cap_at_northpole_grids_ddk = MakeGridDH(filt_SHC_ddk,sampling=2)
        cap_at_northpole_grids_gau = MakeGridDH(filt_SHC_gau,sampling=2)
        
        lons,lats = cap_at_northpole_grids_class.lons(),cap_at_northpole_grids_class.lats()
        fig_name = 'raw_ddk_gau.png'
        fig_name1,fig_name2,fig_name3 = 'psf_raw.png','psf_ddk.png','psf_gau.png'
        magnify = 1e3
        plot_at_northpole(lons,lats,cap_at_northpole_grids,fig_dir+fig_name1,'[mm w.e.]',magnify)
        plot_at_northpole(lons,lats,cap_at_northpole_grids_ddk,fig_dir+fig_name2,'[mm w.e.]',magnify)
        plot_at_northpole(lons,lats,cap_at_northpole_grids_gau,fig_dir+fig_name3,'[mm w.e.]',magnify) 

        fig = plt.figure(dpi=200)
        ax = fig.add_subplot(111)
    
        ax.plot(thetas,value*magnify,'y',label='raw')
        ax.plot(thetas, value_ddk*magnify, 'b--',label='ddk')
        ax.plot(thetas, value_gau*magnify, 'r--',label='gaussian')
        ax.set_xlabel(r'$\theta$ [deg]',size = 'small')
        ax.set_ylabel('[mm w.e.]',size = 'small')
        ax.legend() 
        plt.savefig(fig_dir+fig_name,bbox_inches='tight')           

        