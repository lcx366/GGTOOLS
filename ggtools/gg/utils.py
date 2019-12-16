import numpy as np
import pySphericalPolygon as pysp
from pyshtools.spectralanalysis import Curve2Mask

from ..ggclasses.class_Nodes import Nodes
        

def print_error_grace(source, D, RL):
    '''
    If source, D, and RL do not meet the input requirements, an error message will be printed.
    
    Usage: 
    print_error(source, D, RL)
    
    Inputs:
    source -> [str] source of the GRACE Level-2 solutions; the available sources are CSR, GFZ, JPL
    D -> [int] degree of the GRACE Level-2 solutions; the degree to choose from is either 60 or 96
    RL -> [str] release of the GRACE Level-2 solutions; currently only RL06 is available

    Outputs: raise error message
        
    Examples:
    >>> print_error_grace('CNES/GRGS', 60, 'RL06')
    Exception: Currently, only data from institutions such as CSR, GFZ, and JPL are feasible
    >>> print_error_grace('CSR', 100, 'RL06')
    ValueError: Degree should be either 60 or 96
    >>> print_error_grace('GFZ', 96, 'RL05')
    Exception: Currently, only release RL06 is feasible
    '''
    if source not in ['CSR','GFZ','JPL']:
        raise Exception('Currently, only data from institutions such as CSR, GFZ, and JPL are feasible')
        
    if D not in [60,96]:
        raise ValueError('Degree should be either 60 or 96')  
        
    if RL != 'RL06':
        raise Exception('Currently, only release RL06 is feasible')
    
    return None

def print_error_gldas(source, res):
    '''
    If source and res do not meet the input requirements, an error message will be printed.
    
    Usage: 
    print_error(source, D, RL)
    
    Inputs:
    source -> [str] source of the GLDAS solutions; available option is 'NOAH'
    res -> [str] space resolution of the GLDAS grid data; available options are '1deg' and '0.25deg'

    Outputs: raise error message
    '''
    if source not in ['NOAH']:
        raise Exception('Currently, only data from NOAH are feasible')
        
    if res == '0.25deg':
        dir_res = '025'
        val_res = float(res.strip('deg'))
    elif res == '1deg':
        dir_res = '10'
        val_res = float(res.strip('deg'))
    else:
        raise Exception('Currently, avaliable resolution are 1deg and 0.25deg')
    
    return val_res,dir_res

def month2int(month_list):
    '''
    Given a list of month list, translate it to an array of month sequence.
    
    Usage: 
    month_seq = month2int(month_list)
    
    Inputs:
    month_list -> [str array] list of month, e.g., ['2012-03','2012-04','2012-06','2012-10']

    Outputs:
    month_seq -> [int array] month sequence

    Examples:
    >>> month_list = ['2012-03','2012-04','2012-06','2012-10']
    >>> month_seq = month2int(month_list)
    >>> print(month_seq)
    [0 1 3 7]
    '''
    month_num = []
    n = len(month_list)
    for i in range(n):
        month_num.append((int(month_list[i][:4]) - int(month_list[0][:4]))*12 + (int(month_list[i][5:]) - int(month_list[0][5:])))
    return np.array(month_num)

def med(x):
    '''
    Calculate the middle divisor. This program is used to determine the optimal grid line spacing for regional maps.

    Usage: 
    y = med(x)
    
    Inputs:
    x -> [int] The longitude or latitude span of a map. It cannot be a prime number.

    Outputs:
    y -> [int] the optimal grid line spacing

    Examples:
    >>> print(med(20))
    4
    >>> print(med(25))
    5
    '''
    y = np.unique(np.gcd(np.arange(x),x))
    n = len(y)
    if n%2 == 1:
        return y[n//2]
    else:
        return y[n//2-1]

def solid_angle_ratio(theta_a,theta_b):
    '''
    Calculate the ratio of a solid angle of an ellipse on a sphere to 4pi.

    Usage: 
    ratio = solid_angle_ratio(theta_a,theta_b)
    
    Inputs:
    theta_a -> [float] Semi-minor axis of ellipse
    theta_b -> [float] Semi-major axis of ellipse

    Outputs:
    ratio -> [float] ratio of the solid angle of the ellipse to 4pi
    '''
    return np.sin(np.deg2rad(theta_a+theta_b)/4)**2         

def crop_region(lons,lats,grids,grids_std,region):
    '''
    Crop the global grid data to a region of interested.

    Usage: 
    lons_region,lats_region,grids_region,grids_std_region,lons_flag,lats_flag = crop_region(lons,lats,grids,grids_std,region)
    
    Inputs:
    lons -> [float array] longitudes of the global grid data
    lats -> [float array] latitudes of the global grid data
    grids -> [float 3d array] global grids
    grids_std -> [float 3d array] standard deviation for global grids
    region -> [int array or list] range of longitude and latitude. It should be defined as [leftlon, rightlon, lowerlat, upperlat]

    Outputs:
    lons_region -> [float array] longitudes of the cropped grid data
    lats_region -> [float array] latitudes of the cropped grid data
    grids_region -> [float 3d array] cropped grids
    grids_std_region -> [float 3d array] standard deviation for cropped grids
    lons_flag -> [bool array] Identify which longitude values are within the area of interest. If Ture, inside. If False, outside.
    lats_flag -> [bool array] Identify which latitude values are within the area of interest. If Ture, inside. If False, outside.
    '''
    leftlon, rightlon, lowerlat, upperlat = region
    buffer_edge = lons[1] - lons[0]

    llon_buf = leftlon - buffer_edge
    rlon_buf = rightlon + buffer_edge
    llat_buf = lowerlat - buffer_edge
    ulat_buf = upperlat + buffer_edge
        
    # longitudes and latitudes within the region
    lons_flag = np.logical_and(lons >= llon_buf, lons <= rlon_buf)
    lats_flag = np.logical_and(lats >= llat_buf, lats <= ulat_buf)

    lons_region = lons[lons_flag]
    lats_region = lats[lats_flag]
        
    boundary_lons, = np.diff(lons_flag).nonzero()
    boundary_lats, = np.diff(lats_flag).nonzero()
    grids_region = grids[:,boundary_lats[0]+1:boundary_lats[-1]+1,boundary_lons[0]+1:boundary_lons[-1]+1]
    grids_std_region = grids_std[:,boundary_lats[0]+1:boundary_lats[-1]+1,boundary_lons[0]+1:boundary_lons[-1]+1]

    return lons_region,lats_region,grids_region,grids_std_region,lons_flag,lats_flag 

def yx2latlon(lats_region,lons_region,yx):
    '''
    Transform index to latitudes and longitudes.
 
    Usage: 
    latlon = yx2latlon(lats_region,lons_region,yx)
 
    Inputs:
    lats_region -> [float array] latitudes of a region
    lons_region -> [float array] longitudes of a region
    yx -> [float 2d arrays] index of points in the form of [[y0,x0],..[yn,xn]]; note that the y coordinates are counted from top to bottom;
 
    Outputs:
    latlon -> [float arrays] latitudes and longitudes in the form of [[lat0,lon0],..[latn,lonn]]

    Example:
    >>> lats_region = np.arange(20,50)
    >>> lons_region = np.arange(70,110)
    >>> yx = np.array([[1,2],[13,17],[24,36]])
    >>> latlon = yx2latlon(lats_region,lons_region,yx)
    >>> print(latlon)
    [[ 21  72]
     [ 33  87]
    [ 44 106]]
    '''
    ys,xs = yx[:,0],yx[:,1]
    step_lats,step_lons = lats_region[1]-lats_region[0],lons_region[1]-lons_region[0]
    lats = lats_region[0]+ys*step_lats
    lons = lons_region[0]+xs*step_lons
    return np.stack((lats,lons),axis=1)

def latlon2yx(lats_region,lons_region,latlon):
    '''
    Transform latitudes and longitudes to index.
 
    Usage:
    yx = latlon2yx(lats_region,lons_region,latlon)
 
    Inputs:
    lats_region -> [float array] latitudes of a region
    lons_region -> [float array] longitudes of a region
    latlon -> [float arrays] latitudes and longitudes in the form of [[lat0,lon0],..[latn,lonn]]
 
    Outputs:
    yx: -> [float 2d arrays] index of points in the form of [[y0,x0],..[yn,xn]]; note that the y coordinates are counted from top to bottom
    '''
    lats,lons = latlon[:,0],latlon[:,1]
    step_lats,step_lons = lats_region[1]-lats_region[0],lons_region[1]-lons_region[0]
    ys = np.round((lats - lats_region[0])/step_lats)
    xs = np.round((lons - lons_region[0])/step_lons)
    
    return np.stack((ys.astype(int),xs.astype(int)),axis=1)    

def generate_nodes(grids,points,research_boundary,mode=None):
    '''
    Create nodes within a study area using a set of points or polygons.
 
    Usage:
    nodes,nodes_index = generate_nodes(lats,lons,points,research_boundary)
    nodes,nodes_index = enerate_nodes(lats,lons,points,research_boundary,mode='polygon')
 
    Inputs:
    grids -> [object] instance of class Grid
    points -> [float 2d/3d arrays] latitudes and longitudes in form of [[lat0,lon0],..[latn,lonn]] or [[[lat0,lon0],..[latn,lonn]],[[lat0,lon0],..[latn,lonn]]]
    If 2d, these points can be either a series of discrete points or a polygon. If 3d, these points make up multiple polygons. 
    research_boundary -> [float 2d array] Study area surrounded by a polygon, which is defined as a series of points [[lat0,lon0],..,[latn,lonn]].

    Parameters:
    mode -> [optional, str, default = 'polygon'] If None, mascons are created by a set of discrete points. If 'polygon', mascons are created by a polygon or multiple polygons.

    Outputs:
    nodes -> [float 2d arrays] grid nodes of mascons in form of [[lat0,lon0],..[latn,lonn]]
    nodes_index -> [float 2d arrays] index of nodes in the form of [[y0,x0],..[yn,xn]]; note that the y coordinates are counted from top to bottom
    '''
    lats,lons = grids.lats,grids.lons
    nlat,nlon = len(lats),len(lons)
    lmax = int(nlat/2-1)
    lats_step = lats[1] - lats[0]
    lons_step = lons[1] - lons[0]
    layer = np.zeros((nlat,nlon)).astype(int)
    
    # construct a polygon region from points
    study_area = pysp.SpericalPolygon(research_boundary)
        
    if mode is None:
        nodes = np.zeros_like(points)

        for i in range(len(points)):
            nodes_lat_flag = np.abs(points[i,0] - lats) < -lats_step/2
            nodes_lon_flag = np.abs(points[i,1] - lons) < lons_step/2
            
            lat_index, = np.where(nodes_lat_flag)
            lon_index, = np.where(nodes_lon_flag)
            node_lat,node_lon = lats[nodes_lat_flag],lons[nodes_lon_flag]
            nodes[i] = [node_lat,node_lon]

        nodes = np.unique(nodes,axis=0)   

        # identify whether these nodes are in the region
        inside_flag = study_area.contains_points(nodes) 
        nodes_inside = nodes[inside_flag]   
        nodes_inside_index = latlon2yx(lats,lons,nodes_inside)
        
        return Nodes(nodes_inside,nodes_inside_index,grids.region) 
          
    elif mode == 'polygon':
        if points.ndim == 3:
            for p in points:
                layer = layer | Curve2Mask(nlat,p, 0, sampling = 2)  
        elif points.ndim == 2:
            layer = Curve2Mask(nlat,points, 0, sampling = 2)
        else:
            raise Exception('dimension of polygon should be 2 and polygons should be 3')       
    else:
        raise Exception('only polygon and polygons are avaliable') 
        
    layer = layer & Curve2Mask(nlat,research_boundary, 0, sampling = 2)    
            
    lat_index,lon_index = np.where(layer)
    
    node_lat,node_lon = lats[lat_index],lons[lon_index]
    nodes = np.stack((node_lat,node_lon),axis=1)
    nodes_index = np.stack((lat_index,lon_index),axis=1)
    
    return Nodes(nodes,nodes_index,grids.region)                   