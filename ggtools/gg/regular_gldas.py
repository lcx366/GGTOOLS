import numpy as np
from scipy.interpolate import RectSphereBivariateSpline

def regular_gldas(lats,lons,grids):
    '''
    Normalize the GLDAS grid data to meet the requirements of spherical harmonic expansion with pyshtools based on the sampling theorem of Driscoll and Healy (1994).
    The normalized grid data has the following characteristics: 
    (1) the first latitudinal band corresponds to 90 N, and the latitudinal band for 90 S is not included, and the latitudinal sampling interval is 180/n degrees. 
    (2) the first longitudinal band is 0 E, and the longitude band for 360 E is not included, and the longitudinal sampling interval is 360/n for an equally spaced grid.

    Usage:
    gldas_new = regular_gldas(lats,lons,grids)

    Inputs:
    lats -> [float array] latitudes of gldas grid data
    lons -> [float array] longitudes of gldas grid data
    grids -> [float 2d array] gldas grids data

    Parameters:
    
    Outputs:
    gldas_new -> float 2d array] normalized grids data
    '''
    # Extend the spatial extent of the GLDAS grid to a grid with a shape of n x 2n
    n = int(180/(lats[1] - lats[0]))
    lats_new = np.linspace(lats.max(),-lats.max(),n)
    lons_new = np.linspace(180-lons.max(),180+lons.max(),2*n)

    nlats,nlons = len(lats),len(lons)
    nlats_new,nlons_new = len(lats_new),len(lons_new)

    grids_new = []

    for grid in grids:
        grid_new = np.zeros((nlats_new,nlons_new))
        # Flip the grid data along latitude
        grid_new[:nlats,:] = np.flip(grid,0)
        # Swap the grid data along the longitude direction to make the longitude range from [-180W ~ 180E] to [0E ~ 360E]
        grid_new[:,:nlats_new],grid_new[:,nlats_new:] = grid_new[:,nlats_new:].copy(),grid_new[:,:nlats_new].copy()
        grids_new.append(grid_new)
    grids_new = np.array(grids_new) 
    grids_new[np.isnan(grids_new)] = 0
    
    trunc = np.abs(grids_new[grids_new != 0]).min()/2
    
    lats_interp = lats_new + lons_new.min()
    lons_interp = lons_new - lons_new.min()
    
    grids_interp = []
 
    # colatitude and longitude before interpolation
    colats_rad = np.deg2rad(90 - lats_new)
    lons_rad = np.deg2rad(lons_new)
    
    # colatitude and longitude after interpolation
    colats_rad_interp = np.deg2rad(90 - lats_interp)
    lons_rad_interp = np.deg2rad(lons_interp)
    lons_rad_interp, colats_rad_interp = np.meshgrid(lons_rad_interp, colats_rad_interp)
    
    # grid data after interpolation
    for grid_new in grids_new:
        lut = RectSphereBivariateSpline(colats_rad,lons_rad,grid_new) 
        grid_interp = lut.ev(colats_rad_interp,lons_rad_interp).reshape(nlats_new,nlons_new)
        grids_interp.append(grid_interp)
    grids_interp = np.array(grids_interp) 
    grids_interp[np.abs(grids_interp) < trunc] = 0

    return lats_interp,lons_interp,grids_interp