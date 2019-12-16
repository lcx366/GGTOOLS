import numpy as np
from os import path,makedirs
from urllib.request import urlretrieve
import h5py

from ..ggclasses.class_Mascon import Mascon

def mascon_download(mode=None):
    '''
    Download GRACE mascon data from https://neptune.gsfc.nasa.gov/uploads/grace/mascons_2.4/; if the file to be downloaded is already included in the download directory, the download is automatically skipped.
    
    Usage:
    mascon_download()
    mascon_download('ICE6G')
    
    Parameters:
    mode -> [optional, str, default = None] If None, the standard solution without correction for GIA; comparable to GRACE Level-2 spherical harmonics (GSM) with post-processing corrections applied.
    If 'ICE6G', the standard solution with ICE6G GIA model removed
    
    Outputs: downloaded GRACE mascon data
        
    Examples:
    >>> mascon_download()
    Downloading the GRACE mascon data ... Finished
    '''
    direc = 'MASCON/'
    server = 'https://neptune.gsfc.nasa.gov/uploads/grace/mascons_2.4/'

    if mode is None:
        h5file = 'GSFC.glb.200301_201607_v02.4.h5'
    elif mode == 'ICE6G':
        h5file = 'GSFC.glb.200301_201607_v02.4-ICE6G.h5' 
    else:
        raise Exception('Currently, only ICE6G is available')    

    url = server + h5file   

    if not path.exists(direc): makedirs(direc)
    if not path.exists(direc+h5file):
        print('Downloading the GRACE mascon data',end=' ... ')
        urlretrieve(url, direc+h5file)
        print('Finished')

    # Download the documentation   
    readme_file = 'GSFC_mascons_HDF5_format_v02.4.pdf' 
    url = server + readme_file
    if not path.exists(direc+readme_file): urlretrieve(url, direc+readme_file)

def read_mascon(mode=None):
    '''
    Read the GRACE mascon files. 
    
    Usage: 
    mascon = read_mascon()
    mascon = read_mascon('ICE6G')

    Parameters:
    mode -> [optional, str, default = None] If None, the standard solution without correction for GIA; comparable to GRACE Level-2 spherical harmonics (GSM) with post-processing corrections applied.
    If 'ICE6G', the standard solution with ICE6G GIA model removed
    
    Outputs:
    mascon: instance of Mascon class
        
    Examples:
    >>> mascon = read_mascon()
    >>> print(mascons.area_mascon)
    [12453.61 12453.61 12453.61 ... 12415.49 12415.49 12415.49]
    For more information, please refer to https://neptune.gsfc.nasa.gov/uploads/grace/mascons_2.4/GSFC_mascons_HDF5_format_v02.4.pdf
    '''
    direc = 'MASCON/'

    if mode is None:
        h5file = 'GSFC.glb.200301_201607_v02.4.h5'
    elif mode == 'ICE6G':
        h5file = 'GSFC.glb.200301_201607_v02.4-ICE6G.h5' 
    else:
        raise Exception('Currently, only ICE6G is available') 

    fin = h5py.File(direc+h5file,'r')

    lons = fin['mascon/lon_center'][:].flatten()
    lats = fin['mascon/lat_center'][:].flatten()
    lons_span = fin['mascon/lon_span'][:].flatten()
    lats_span = fin['mascon/lat_span'][:].flatten()
    area_mascon = fin['mascon/area_km2'][:].flatten()
    mascon_cmwe = fin['solution/cmwe'][:]
    mascon_mmwe = mascon_cmwe*10 # convert cm w.e. to mm w.e.

    # close the file
    fin.close()

    return Mascon(lons,lats,lons_span,lats_span,area_mascon,mascon_mmwe)