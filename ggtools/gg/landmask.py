import numpy as np
import xarray as xr
from os import getenv,path,makedirs
from urllib.request import urlretrieve
from pyshtools.shclasses import SHCoeffs,SHGrid

def landmask(lmax,trunc_lat=None):
    '''
    Establish a land window function based on the global terrain data ETOPO5. The land area has a value of 1 and the sea area has a value of 0.
    To set the value of the Antarctica region to 0, just set trunc_lat to -60.

    Usage:
    land_win = land_mask(lmax)
    land_win = land_mask(lmax,-60)

    Inputs:
    lmax -> [int] Degree of the spherical harmonic expansion

    Parameters:
    trunc_lat -> [optional, float, default = None] Truncated latitude. 
    If None, the global land area will be assigned a value of 1, otherwise, only the land area north of the truncated latitude will be assigned a value of 1, and the remaining regions will be assigned a value of 0.
    
    Outputs:
    land_win [bool, 2d array] land window function. The grid satisfies the sampling theorem of Driscoll and Healy (1994).

    For more information, please refer to https://shtools.oca.eu/shtools/public/pyshexpanddh.html
    '''
    # Download ETOPO5 data
    home = getenv('HOME')
    direc = home + '/src/etopo-data/'
    etopo_file = 'etopo5.nc'
    url = 'https://raw.githubusercontent.com/dcherian/tools/master/ROMS/arango/bathymetry/etopo5.nc'
    
    if not path.exists(direc): 
        makedirs(direc)
        print('Downloading the ETOPO5 Earth Surface Topography Data Set',end=' ... ')
        urlretrieve(url, direc + etopo_file)
        print('Finished')

    ds = xr.open_dataset(direc + etopo_file)
    
    # Flip the data along the latitude direction to change the latitude range from [90S ∼ 90N] to [90N ∼ 90S]
    topo_grids = np.flip(np.array(ds['topo']),0)[:-1,:]
    lats = np.flip(np.array(ds['topo_lat']))[:-1]
 
    # Expand grid data to spherical harmonic coefficients
    topo_grids_class = SHGrid.from_array(topo_grids)
    topo_coeffs_class = topo_grids_class.expand(lmax_calc=lmax)
    topo_grids_class = topo_coeffs_class.expand()
    topo_grids = topo_grids_class.data
    
    lats = topo_grids_class.lats()
 
    land_mask = topo_grids > 0 # mask the ocean
    if trunc_lat is not None:
        land_mask[lats < trunc_lat] = False # mask the area south of trunc_lat, i.e., Antarctica
    return land_mask
