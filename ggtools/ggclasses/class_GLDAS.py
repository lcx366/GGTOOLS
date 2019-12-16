import numpy as np
from pyshtools.shclasses import SHGrid
from ..gg.lsm import lsm

from ..gg.utils import crop_region

from .class_Grid import Grid
from .class_GSM import GSM

class GLDAS(object):
    '''
    class GLDAS
    - attributes:
        - info
        - degree_order
        - max_degree
        - max_order
        - title
        - summary
        - resolution
        - institution
        - time_coverage_start
        - time_coverage_end
        - total_month
        - total_month_counts
        - solution_month
        - solution_counts
        - missing_month
        - missing_month_counts
        - missing_solution_flag
        - missing_value
        - date_issued
        - tavg
        - acc
        - inst 
        - filter
        - lons
        - lats
        - data
    
    - methods:
        - twsc_grid
        - twsc_shc
    '''
    def __init__(self,info,lons,lats,data):
        
        self.info = info
        self.degree_order = info['degree_order']
        self.max_degree = info['max_degree']
        self.max_order = info['max_order'] 
        self.title = info['title']
        self.summary = info['summary']
        self.resolution = info['resolution']
        self.institution = info['institution']
        self.time_coverage_start = info['time_coverage_start']
        self.time_coverage_end = info['time_coverage_end']
        self.total_month = info['total_month']
        self.total_month_counts = info['total_month_counts']
        self.solution_month = info['solution_month'] 
        self.solution_counts = info['solution_counts'] 
        self.missing_month = info['missing_month']
        self.missing_month_counts = info['missing_month_counts']
        self.missing_solution_flag = info['missing_solution_flag'] 
        self.missing_value = info['missing_value']
        self.tavg = info['tavg']
        self.acc = info['acc']
        self.inst = info['inst']
        self.date_issued = info['date_issued']
        self.filter = info['filter']
        
        self.lons,self.lats = lons,lats
        self.data = data
        
    def __repr__(self):
    
        return 'title = {:s}\nresolution = {:s}\ninstitution = {:s}\ntime_coverage_start = {:s}\ntime_coverage_end = {:s}\nsolution_counts = {:d}\ntotal_month_counts = {:d}\nmissing_month_counts = {:d}'.format\
        (self.title,self.resolution,self.institution,self.time_coverage_start,self.time_coverage_end,self.solution_counts,self.total_month_counts,self.missing_month_counts)

    def twsc_grid(self,region,mode='classic'):

        info,lats,lons,twscs = lsm(self.lons,self.lats,self.data,'twsc',mode)
        # 计算陆地水储量变化
        info['region'] = region
        info = {**self.info.copy(),**info}
        twscs_std = np.zeros_like(twscs) 
        lons_region,lats_region,twscs_region,twscs_std_region,lons_flag,lats_flag = crop_region(lons,lats,twscs,twscs_std,region)

        return Grid(info,twscs_region,twscs_std_region,lons_region,lats_region,lons_flag,lats_flag)      
    
    def twsc_shc(self,lmax,mode='classic'):

        info,lats,lons,twscs = lsm(self.lons,self.lats,self.data,'twsc',mode)
        # 计算陆地水储量变化
        info = {**self.info.copy(),**info}

        twscs_coeffs = []
        degree_order = self.degree_order

        if degree_order < lmax:
            for twsc in twscs:
                twsc_grids_class = SHGrid.from_array(twsc)
                twsc_coeffs_class = twsc_grids_class.expand()
                twsc_coeffs = np.zeros((2,lmax+1,lmax+1))
                twsc_coeffs[:,:degree_order+1,:degree_order+1] = twsc_coeffs_class.coeffs
                twscs_coeffs.append(twsc_coeffs)

        else:
            for twsc in twscs:
                twsc_grids_class = SHGrid.from_array(twsc)
                twsc_coeffs_class = twsc_grids_class.expand(lmax_calc=lmax)
                twsc_coeffs = twsc_coeffs_class.coeffs
                twscs_coeffs.append(twsc_coeffs)       

        twscs_coeffs = np.array(twscs_coeffs) 
        twscs_coeffs_std = np.zeros_like(twscs_coeffs)

        info['degree_order'] = lmax

        return GSM(info,twscs_coeffs,twscs_coeffs_std)   