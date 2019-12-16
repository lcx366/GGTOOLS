import numpy as np
from numpy.linalg import norm

from pyshtools import SHCoeffs,SHGrid,Slepian
from pyshtools.expand import SHExpandDH,MakeGridDH
from pyshtools.spectralanalysis import SlepianCoeffsToSH,Curve2Mask
from pyshtools.shio import SHCilmToVector

from .filter import filter_gaussian
from .lcurve import L_curve

def scale_factor(shcs,research_boundary,r):
    '''
    Estimate the scale factor(gain factor) used in signal leakage correction. 

    Usage:
    k = scale_factor(SHCs,research_boundary,r)

    Inputs:
    SHCs -> [float 4d array] Spherical harmonic coefficients for TWSC from GLDAS
    research_boundary -> [float 2d array] Study area surrounded by a polygon, which is defined as a series of points [[lat0,lon0],..,[latn,lonn]].
    r -> [float] Gaussian filter radius

    Parameters:
    
    Outputs:
    k -> [float] Scale factor or Gain factor
    '''
    shcs_filter = shcs.filter_gaussian(r)
    series = shcs.study_area(research_boundary)
    series_filter = shcs_filter.study_area(research_boundary)
    k = np.dot(series_filter.qs,series.qs)/norm(series_filter.qs)**2
    return k

def forward_model_initial(FM0_grid,r):
    '''
    Expand the grid data to spherical harmonic coefficients and perform Gaussian filtering, then transfer it back to the grid data. 

    Usage:
    FM0_grid_gau = forward_model_initial(FM0_grid,r)

    Inputs:
    FM0_grid -> [float 2d array] initial grid data
    r -> [float] Gaussian filter radius
    
    Outputs:
    FM0_grid_gau -> [float 2d array] grid data after Gaussian filtering
    '''
    FM0_shc = SHExpandDH(FM0_grid,sampling=2)
    # set deg 0 and deg 1 to zero
    FM0_shc[:,:2,:] = 0
    FM0_shc_gau = filter_gaussian(r,FM0_shc)
    FM0_grid_gau = MakeGridDH(FM0_shc_gau,sampling=2)
    return FM0_grid_gau

def forward_model(FM0_grids,r,amplify_factor = 1.2):
    '''
    Iterative Forward Modeling used to perform signal leakage correction in GRACE data processing.
    
    Usage:
    FM0_grids_iter,FM0_grids_std_iter = forward_model(FM0_grids,r,)
    FM0_grids_iter,FM0_grids_std_iter = forward_model(FM0_grids,r,1.5)

    Inputs:
    FM0_grid -> [float 2d array] initial grid data
    r -> [float] Gaussian filter radius

    Parameters:
    amplify_factor -> [optional, float, default = 1.2] Magnification factor used to accelerate iterative convergence. 
    The larger the value, the faster the convergence, but the lower the accuracy.
    
    Outputs:
    FM0_grids_iter -> [float 2d array] grid data after signal leakage correction
    FM0_grids_std_iter -> [float 2d array] standard deviation of grid data after signal leakage correction. 
    '''
    FM0_grids_iter = []
    
    for FM0_grid in FM0_grids: 
        FM0_grid_iter = FM0_grid
        FM1_grid = forward_model_initial(FM0_grid_iter,r) 
        
        while np.abs(FM0_grid-FM1_grid).max() > 10:
            FM0_grid_iter = FM0_grid_iter + (FM0_grid-FM1_grid)*amplify_factor
            FM1_grid = forward_model_initial(FM0_grid_iter,r)
        FM0_grids_iter.append(FM0_grid_iter)
    return FM0_grids_iter,np.zeros_like(FM0_grids_iter)

def space_domain(grids,nodes,research_boundary,r):  
    '''
    Space domain method used to perform signal leakage correction in GRACE data processing.
    
    Usage:
    Ms,Ms_std = space_domain(grids,nodes,research_boundary,r)

    Inputs:
    grids -> [object] instance of class Grid
    nodes -> [object] grid nodes that represent mascons in the study area. The nodes are described as a set of discerte grid points [[lat0,lon0],..,[latn,lonn]].
    research_boundary -> [float 2d array] Study area surrounded by a polygon, which is defined as a series of points [[lat0,lon0],..,[latn,lonn]].
    r -> [float] Gaussian filter radius

    Outputs:
    Ms -> [float 2d array] series of mascons or rate of mascons. The first dimension is time and the second dimension is the value of mascon. 
    If the size of the first dimension is 1, it means rate of mascon.
    Ms_std -> [float 2d array] standard deviation for series of mascons or rate of mascons.
    '''
    
    psf_grids_gau,mas = [],[]
    lats_flag,lons_flag = grids.lats_flag,grids.lons_flag
    nlat = len(lats_flag)
    lmax = int(nlat/2-1)

    if grids.region != nodes.region:
        raise Exception('The range of the nodes should be consistent with the range of the grid.')

    for node in nodes.nodes:

        psf_grid_class = SHGrid.from_cap(0.1,node[0],node[1],lmax)
        psf_grid =psf_grid_class.data

        psf_shc = SHExpandDH(psf_grid,sampling=2) 
        psf_shc_gau = filter_gaussian(r,psf_shc)
        psf_grid_gau = MakeGridDH(psf_shc_gau,sampling=2)
    
        psf_grids_gau.append(psf_grid_gau)  
    psf_grids_gau = np.array(psf_grids_gau)  

    mask_boundary = Curve2Mask(nlat, research_boundary, 0, sampling = 2)
    psf_grids_simply = psf_grids_gau[:,mask_boundary.astype(bool)]
    
    if grids.region != 'globe':
        mask_boundary = mask_boundary[lats_flag][:,lons_flag]
        
    grids_simply = grids.grids[:,mask_boundary.astype(bool)]
        
    A = psf_grids_simply.T
    for grid_simply in grids_simply:
        lamb,ma = L_curve(A,grid_simply)
        mas.append(ma)
    mas = np.array(mas)    
    return mas,np.zeros_like(mas) 

def spectral_domain(shcs,nodes,research_boundary,r,mode=None,ratio=None):
    '''
    Spectrum domain method used to perform signal leakage correction in GRACE data processing.
    
    Usage:
    Ms,Ms_std = spetral_domain(SHCs,nodes,research_boundary,r)
    Ms,Ms_std = spetral_domain(SHCs,nodes,research_boundary,r,'window',ratio)

    Inputs:
    shcs -> [float 4d array] spherical harmonic coefficients before signal leakage correction
    nodes -> [object] grid nodes that represent mascons in the study area. The nodes are described as a set of discerte grid points [[lat0,lon0],..,[latn,lonn]].
    research_boundary -> [float 2d array] Study area surrounded by a polygon, which is defined as a series of points [[lat0,lon0],..,[latn,lonn]].
    r -> [float] Gaussian filter radius

    Parameters:
    mode -> [optional, str, default = None] If None, the global spherical harmonic coefficients are used to fit mascons. If 'window', the windowed spherical harmonic coefficients are used to fit mascons.
    ratio -> [optional, str, default = None] The ratio of the study area to the global area, which is mainly used to determine the Shannon number approximately. If None, the mode must be set to None. 

    Outputs:
    Ms -> [float 2d array] series of mascons or rate of mascons. The first dimension is time and the second dimension is the value of mascon. 
    If the size of the first dimension is 1, it means rate of mascon.
    Ms_std -> [float 2d array] standard deviation for series of mascons or rate of mascons. 
    '''
    lmax = shcs.shape[-1]-1
    nlat = (lmax+1)*2
    psf_shcs_gau,mas = [],[]
    
    for node in nodes.nodes:
        
        psf_grid_class = SHGrid.from_cap(0.1,node[0],node[1],lmax)
        psf_grid =psf_grid_class.data
        
        psf_shc = SHExpandDH(psf_grid,sampling=2) 
        psf_shc_gau = filter_gaussian(r,psf_shc)

        psf_shcs_gau.append(SHCilmToVector(psf_shc_gau))  
    psf_shcs_gau = np.array(psf_shcs_gau) 

    A = psf_shcs_gau.T

    if mode is None:
        for shc in shcs:
            y = SHCilmToVector(shc) 
            lamb,ma = L_curve(A,y)
            mas.append(ma)
    else:
        mask_boundary = Curve2Mask(nlat, research_boundary, 0, sampling = 2)

        # estimate the Shannon number
        N = np.ceil((lmax+1)**2*ratio) # Rough Shannon number
        slepian = Slepian.from_mask(mask_boundary,lmax,N)
        slepian_tapers = slepian.tapers
        N = np.ceil(slepian.shannon) # Accurate Shannon Number
    
        for shc in shcs:
            shc_class = SHCoeffs.from_array(shc)
            slepian_coeffs = slepian.expand(shc_class,N)
            shc_slepian = SlepianCoeffsToSH(slepian_coeffs.falpha,slepian_tapers,N)
    
            y = SHCilmToVector(shc_slepian) 
            lamb,ma = L_curve(A,y)
            mas.append(ma)
    mas = np.array(mas)    
    
    return mas,np.zeros_like(mas)     