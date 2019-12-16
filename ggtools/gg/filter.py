import numpy as np
from oct2py import octave
import pkg_resources

def filter_ddk(filter_type,shc,shc_std=None):
    '''
    DDK filter used to attenuate noise described as striping patterns in GRACE GSM data. 
    According to the filtering strength, there are 8 kinds of ddk filter. 
    From DDK1 to DDK8, the filtering effect gradually weakens.

    Usage:
    filt_SHC = filt_DDK('DDK3',SHC)
    filt_SHC,filt_SHC_std = filt_DDK('DDK4',SHC,SHC_std)

    Inputs:
    filt_type -> [str] Types of DDK filter. Available options are 'DDK1' to 'DDK8'.
    SHC -> [float 3d/4d array] Fully normalized spherical harmonics coefficients(Stokes coefficients). The dimension of SHC can either be 3d or 4d.
    If 3d, its shape is (i,l,m), where i = 0 for Clm and i = 1 for Slm. If 4d, its shape is (k,i,l,m).

    Parameters:
    SHC_std -> [optional, float 3d/4d array, default = None] Standard deviation for SHC. 
    If None, the standard deviation for the filtered SHC will not be estimated, and only the filtered SHC is returned.

    Outputs:
    filt_SHC -> [float 3d/4d array] filtered SHC
    filt_SHC_std -> [float 3d/4d array] standard deviation for filtered SHC

    For more information, please refer to https://github.com/strawpants/GRACE-filter
    '''  
    filter_shc = np.zeros_like(shc)
    filter_shc_std = np.zeros_like(shc_std)
        
    # add the matlab scripts path
    scripts_path = pkg_resources.resource_filename('ggtools', 'grace_filter/src/matlab') 
    octave.addpath(scripts_path)   
    # read the filter matrix
    filter_path = pkg_resources.resource_filename('ggtools', 'grace_filter/data/DDK/')
    if filter_type == 'DDK1':
        Wbd = octave.read_BIN(filter_path + 'Wbd_2-120.a_1d14p_4')
    elif filter_type == 'DDK2':
        Wbd = octave.read_BIN(filter_path + 'Wbd_2-120.a_1d13p_4')
    elif filter_type == 'DDK3':
        Wbd = octave.read_BIN(filter_path + 'Wbd_2-120.a_1d12p_4')
    elif filter_type == 'DDK4':
        Wbd = octave.read_BIN(filter_path + 'Wbd_2-120.a_5d11p_4')
    elif filter_type == 'DDK5':
        Wbd = octave.read_BIN(filter_path + 'Wbd_2-120.a_1d11p_4')
    elif filter_type == 'DDK6':
        Wbd = octave.read_BIN(filter_path + 'Wbd_2-120.a_5d10p_4')
    elif filter_type == 'DDK7':
        Wbd = octave.read_BIN(filter_path + 'Wbd_2-120.a_1d10p_4')
    elif filter_type == 'DDK8':
        Wbd = octave.read_BIN(filter_path + 'Wbd_2-120.a_5d9p_4')
    else:
        raise Exception('Currently, only DDK1~DDK8 are feasible.')
            
    if shc.ndim == 4:
            
        if shc_std is None:
            for i in range(shc.shape[0]):
                filter_shc[i,0],filter_shc[i,1] = octave.filterSH(Wbd,shc[i,0],shc[i,1],nout=2) 
            return filter_shc 
        else:
            for i in range(shc.shape[0]):
                filter_shc[i,0],filter_shc[i,1],filter_shc_std[i,0],filter_shc_std[i,1] = octave.filterSH(Wbd,shc[i,0],shc[i,1],shc_std[i,0],shc_std[i,1],nout=4) 
            return filter_shc,filter_shc_std
            
    elif shc.ndim == 3: 
            
        if shc_std is None:
            filter_shc[0],filter_shc[1] = octave.filterSH(Wbd,shc[0],shc[1],nout=2) 
            return filter_shc
        else:
            filter_shc[0],filter_shc[1],filter_shc_std[0],filter_shc_std[1] = octave.filterSH(Wbd,shc[0],shc[1],shc_std[0],shc_std[1],nout=4) 
            return filter_shc,filter_shc_std
    else:
        raise Exception('Dimension of the SHC data is not correct. It should be 3 or 4.')

def filter_gaussian(r,shc,shc_std = None,a = 6378.1363):
    '''
    Gaussian filter used to attenuate noise described as striping patterns in GRACE GSM data. 
    The effect of filtering can be changed by adjusting the filtering radius.
    Larger filter radius can significantly suppress noise, but at the same time attenuate the signals.

    Usage:
    filt_SHC = filt_gaussian(SHC) 
    filt_SHC = filt_gaussian(SHC,r = 300)
    filt_SHC,filt_SHC_std = filt_gaussian(SHC,SHC_std)
    filt_SHC,filt_SHC_std = filt_gaussian(SHC,SHC_std,r = 400)

    Inputs:
    r -> [float] Gaussian filer radius in km
    SHC -> [float 3d/4d array] Fully normalized spherical harmonics coefficients(Stokes coefficients). The dimension of SHC can either be 3d or 4d.
    If 3d, its shape is (i,l,m), where i = 0 for Clm and i = 1 for Slm. If 4d, its shape is (k,i,l,m).

    Parameters:
    SHC_std -> [optional, float 3d/4d array, default = None] Standard deviation for SHC. 
    If None, the standard deviation for the filtered SHC will not be estimated, and only the filtered SHC is returned.
    a -> [optional, float, default = 6378.1363] Average Earth radius in km

    Outputs:
    filt_SHC -> [float 3d/4d array] filtered SHC
    filt_SHC_std -> [float 3d/4d array] standard deviation for filtered SHC
    '''   
    n = shc.shape[-1]-1
    b = np.log(2)/(1-np.cos(r/a))
    W = np.ones(n+1)
    W[1] = (1+np.exp(-2*b))/(1-np.exp(-2*b))-1/b
        
    for i in range(1,n):
        W[i+1] = -(2*i+1)/b*W[i]+W[i-1]

    filter_shc = shc*W[:,None]
        
    if shc_std is None:
        return filter_shc
    else:
        filter_shc_std = shc_std*W[:,None]
        return filter_shc,filter_shc_std           

def filter_gaussian_inverse(r,shc,shc_std=None,a = 6378.1363):
    '''
    Inversion of Gaussian filter. It is used to recover the signal in the process of leakage correction in GRACE GSM data. 

    Usage:
    recover_SHC = filt_gaussian_inverse(SHC) 
    recover_SHC = filt_gaussian_inverse(SHC,r = 300)
    recover_SHC,recover_SHC_std = filt_gaussian_inverse(SHC,SHC_std)
    recover_SHC,recover_SHC_std = filt_gaussian_inverse(SHC,SHC_std,r = 400)

    Inputs:
    r -> [float] Gaussian filer radius in km
    SHC -> [float 3d/4d array] Fully normalized spherical harmonics coefficients(Stokes coefficients). The dimension of SHC can either be 3d or 4d.
    If 3d, its shape is (i,l,m), where i = 0 for Clm and i = 1 for Slm. If 4d, its shape is (k,i,l,m).

    Parameters:
    SHC_std -> [optional, float 3d/4d array, default = None] Standard deviation for SHC. 
    If None, the standard deviation for the recovered SHC will not be estimated, and only the recovered SHC is returned.
    a -> [optional, float, default = 6378.1363] Average Earth radius in km

    Outputs:
    recover_SHC -> [float 3d/4d array] recovered SHC
    recover_SHC_std -> [float 3d/4d array] standard deviation for recovered SHC
    ''' 
    n = shc.shape[-1]-1
    b = np.log(2)/(1-np.cos(r/a))
    W = np.ones(n+1)
    W[1] = (1+np.exp(-2*b))/(1-np.exp(-2*b))-1/b
        
    for i in range(1,n):
        W[i+1] = -(2*i+1)/b*W[i]+W[i-1]

    recover_shc = shc/W[:,None]
        
    if shc_std is None:
        return recover_shc
    else:
        recover_shc_std = shc_std/W[:,None]
        return recover_shc,recover_shc_std            