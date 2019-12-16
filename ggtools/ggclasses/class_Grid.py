from os import path,makedirs
import numpy as np

from pyshtools.spectralanalysis import Curve2Mask

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from ..gg.utils import month2int,med,crop_region
from ..gg.lsq import lsqm,ilsqm,wlsqm,iwlsqm
from ..gg.leakage import forward_model,space_domain

from .class_Series import Series

class Grid(object):
    '''
    class Grid
    - attributes
        - info -> All information about the object
        - degree_order 
        - max_degree
        - max_order
        - normalization
        - permanent_tide
        - earth_gravity_param
        - mean_equator_radius
        - background_gravity
        - title
        - summary
        - institution
        - processing_level
        - product_version
        - time_coverage_start
        - time_coverage_end
        - total_month -> Months over a time interval regardless of the existence of the solutions
        - total_month_counts
        - solution_month
        - solution_counts
        - missing_month
        - missing_month_counts
        - missing_solution_flag -> If True, the monthly solution is missing, otherwise, the monthly solution exists
        - unused_days -> Unused days for monthly solutions
        - date_issued
        - equi_material -> Equivalent material used to represent mass per unit area
        - filter -> filter applied to monthly solutions
        - region
        - lons
        - lats
        - lons_flag
        - lats_flag
        - grids
        - grids_std
    
    - methods:
        - rate
        - study_area
        - plot
        - set region
        - leakage_correction
    '''
    def __init__(self,info,grids,grids_std,lons,lats,lons_flag,lats_flag):
        
        self.info = info
        self.degree_order = info['degree_order']
        self.max_degree = info['max_degree']
        self.max_order = info['max_order'] 
        self.normalization = info['normalization'] 
        self.permanent_tide = info['permanent_tide'] 
        self.earth_gravity_param = info['earth_gravity_param']
        self.mean_equator_radius = info['mean_equator_radius']
        self.background_gravity = info['background_gravity']
        self.title = info['title']
        self.summary = info['summary']
        self.institution = info['institution']
        self.processing_level = info['processing_level']
        self.product_version = info['product_version']
        self.time_coverage_start = info['time_coverage_start']
        self.time_coverage_end = info['time_coverage_end']
        self.total_month = info['total_month']
        self.total_month_counts = info['total_month_counts'] 
        self.solution_counts = info['solution_counts'] 
        self.solution_month = info['solution_month'] 
        self.missing_month = info['missing_month']
        self.missing_month_counts = info['missing_month_counts']
        self.missing_solution_flag = info['missing_solution_flag']
        self.unused_days = info['unused_days']
        self.date_issued = info['date_issued']
        self.equi_material = info['equi_material']
        self.filter = info['filter']
        self.region = info['region']
        
        self.grids,self.grids_std = grids,grids_std
        self.lons,self.lats = lons,lats
        self.lons_flag,self.lats_flag = lons_flag,lats_flag
        
    def __repr__(self):
    
        return 'title = {:s}\ntime_coverage_start = {:s}\ntime_coverage_end = {:s}\nsolution_counts = {:d}\ntotal_month_counts = {:d}\nmissing_month_counts = {:d}\nequi_material = {:s}\nregion = {:s}'.format\
        (self.title,self.time_coverage_start,self.time_coverage_end,self.solution_counts,self.total_month_counts,self.missing_month_counts,self.equi_material,str(self.region))
    
    def rate(self,mode='ILSQM'):
        '''
        Estimate the annual change rate of grids for Surface Mass Anomaly in Equivalent Water(or Ice, Sand) Thickness(EWT) using the linearly fitting method.
        There are four methods for linearly fitting, including 'LSM', 'ILSM', 'WLSM', and 'IWSLM'. The ILSM is default and recommended.
        
        Usage: xxx_sma_grid_rate = xxx_sma_grid.rate() or xxx_sma_grid_rate = xxx_sma_grid.rate('IWLSM')

        Parameters:
        -----------
            lsm [str] [optional, default: ILSM] alternertively, 'LSM', 'ILSM', 'WLSM', and 'IWSLM' are available, where
            'LSM'   --  Least Square Method
            'ILSM'  --  Iterative Least Square Method
            'WLSM'  --  Weighted Least Square Method
            'IWSLM' --  Iterative Weighted Least Square Method
            
        Returns:
        -----------
            xxx_sma_grid_rate: instance of GRID class
            
        Examples:
        -----------
        >>> csr_gsm = read_GSM('CSR',96)
        >>> gfz_gsm = read_GSM('GFZ',96)
        >>> jpl_gsm = read_GSM('JPL',96)
        >>> slr_c20 = read_SLR_C20()
        >>> comb_gsm = GSM_average([csr_gsm.deaverage(),gfz_gsm.deaverage(),jpl_gsm.deaverage()])
        >>> comb_gsm_r = comb_gsm.replace_slr_c20(slr_c20.deaverage())
        >>> comb_gsm_ddk5 = comb_gsm_r.filt_DDK('DDK5')
        >>> comb_gsm_ddk5_rate = comb_gsm_ddk5.rate()
        >>> print(comb_gsm_ddk5_rate.title)
        'Annual change rate of DDK5 filtered Combined Deaveraged GRACE & GRACE-FO Geopotential Coefficients CSR RL06, GFZ RL06, JPL RL06 with C20 replaced by the SLR measurements'
        >>> print(comb_gsm_ddk5_rate.summary)
        'Spherical harmonic coefficients representing an estimate of annual change rate of the mean gravity field of Earth during the specified timespan derived from GRACE & GRACE-FO mission measurements. These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes. Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product. Note that the 2nd degree terms have been replaced by the C20 values from SLR. Also note that C20 values from SLR also experienced the DDK5 filtering.'
        >>> sma_rate = comb_gsm_ddk5_rate.sma()
        >>> print(sma_rate.title)
        'Stokes coefficients for annual change rate of Surface Mass Anomaly(SMA) in Equivalent Water Thickness(EWT) derived from the Annual change rate of DDK5 filtered Combined Deaveraged GRACE & GRACE-FO Geopotential Coefficients CSR RL06, GFZ RL06, JPL RL06 with C20 replaced by the SLR measurements'
        >>> print(sma_rate.summary)
        'Spherical harmonic coefficients representing an estimate of annual change rate of the Surface Mass Anomaly(SMA) expressed in terms of Equivalent Water[1000kg/m3] Thickness(EWT) with unit of [mm w.e./yr] during the specified timespan derived from GRACE & GRACE-FO mission measurements. These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes. Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product. Note that the 2nd degree terms have been replaced by the C20 values from SLR. Also note that C20 values from SLR also experienced the DDK5 filtering.'
        '''
        grids,grids_std = self.grids,self.grids_std
        lats,lons = self.lats,self.lons
        month = month2int(self.solution_month)
        grids_rate,grids_rate_std = [np.zeros_like(grids[0]) for k in range(2)]
        
        for i in range(len(lats)):    
            for j in range(len(lons)): 
                if mode is 'LSQM':
                    grids_rate[i,j],grids_rate_std[i,j],_,_ = lsqm(month,grids[:,i,j])
                elif mode is 'WLSQM':
                    grids_rate[i,j],grids_rate_std[i,j],_,_ = wlsqm(month,grids[:,i,j],grids_std[:,i,j])
                elif mode is 'ILSQM':
                    grids_rate[i,j],grids_rate_std[i,j],normal,_,_ = ilsqm(month,grids[:,i,j])
                elif mode is 'IWLSQM':
                    grids_rate[i,j],grids_rate_std[i,j],normal,_,_ = iwlsqm(month,grids[:,i,j],grids_std[:,i,j])
                else:
                    raise Exception('Currently, the least square method can only be LSQM, WLSQM, ILSQM, and IWLSQM.')
        info = self.info.copy()
        info['title'] = 'Annual change rate of ' + info['title']
        info['summary'] = info['summary'].replace('an estimate of','an estimate of annual change rate of')
        for em in ['w','i','s']:
            info['summary'] = info['summary'].replace('[mm '+em+'.e.]','[mm '+em+'.e./yr]')

        return Grid(info,np.array([grids_rate]),np.array([grids_rate_std]),lons,lats,self.lons_flag,self.lats_flag)
    
    def study_area(self,points,north_pole=False,central_meridian=False):
        
        lats = self.lats
        lons_flag,lats_flag = self.lons_flag,self.lats_flag
        a = float(self.mean_equator_radius.partition('m')[0])/1e3 # km
        
        if self.equi_material is 'Water':
            rho = 1000
        elif self.equi_material is 'Ice':
            rho = 917
        elif self.equi_material is 'Sand':
            rho = 1442
        
        lat_step = lats[1] - lats[0]
        qs,qs_std = [],[]

        mask = Curve2Mask(2*(self.degree_order+1),points,north_pole,sampling=2,centralmeridian=central_meridian)
        mask_region = mask[lats_flag,:][:,lons_flag]
        
        grid_area = (np.deg2rad(lat_step)*a)**2 # km^2
        grid_metric = np.cos(np.deg2rad(lats))[:,None] # cos metric
        
        area = np.sum(mask_region*grid_metric)*grid_area # km^2
        
        for grid in self.grids:
            q = np.sum(grid*mask_region*grid_metric)*grid_area*rho/1e9 # Gt
            qs.append(q)
        for grid_std in self.grids_std:    
            q_std = np.sqrt(np.sum(grid_std**2*mask_region*grid_metric**2))*grid_area*rho/1e9 # Gt
            qs_std.append(q_std)
        qs,qs_std = np.array(qs),np.array(qs_std) 
        
        info = self.info.copy()
        info['title'] = info['title'].replace('Grids','Sum(over the study area) of grids')
        
        return Series(info,area,qs,qs_std)
       
        
    def plot(self,fig_name,ylabel,block=None,polygons=None,circles=None):
        
        '''
        Plot grids data.

        Usage: plot_figure(region,buffer_edge,lons_region,lats_region,data_region,fig_name,circles)

        INPUTS:
            region: [float array] range of region, for example, [96.0,120.0,21.0,39.0] means 
            the boundaries;
            fig_name: [str] filename of the output figure;
            ylabel: [str] the label of y axis;
            circles: [float array] the boundary of circles or polygons; if None, no circles
            will be plotted.

        OUTPUTS:
            figure: the output figure;
        '''
        
        if not 'rate' in self.title:
            raise Exception('The shape of the grid data to be plotted should be like (1,d1,d2)')
        
        lons_region,lats_region = self.lons,self.lats
        grids_region = self.grids
        img_extent = self.region
        
        leftlon, rightlon, lowerlat, upperlat = img_extent
        lon_step = rightlon - leftlon
        lat_step = upperlat - lowerlat
        
        buffer_edge = lons_region[1] - lons_region[0]
        gridlons = np.append(lons_region-buffer_edge/2,lons_region[-1]+buffer_edge/2)
        gridlats = np.append(lats_region+buffer_edge/2,lats_region[-1]-buffer_edge/2)
        
        fig_dir = 'figures/'
        if not path.exists(fig_dir): makedirs(fig_dir) 
        
        # set the projection to PlateCarree
        proj = ccrs.PlateCarree(central_longitude=180)
        
        # plot
        for grid_region in grids_region:
            plt.clf()
            fig = plt.figure(dpi=200)
            ax = fig.add_subplot(1, 1, 1,projection = proj)

            # set the range of the area of interest
            ax.set_extent(img_extent,crs = ccrs.PlateCarree())
            gl = ax.gridlines(xlocs=np.linspace(leftlon,rightlon,med(int(lon_step))+1),ylocs=np.linspace(lowerlat,upperlat,med(int(lat_step))+1),draw_labels=True,linestyle='--',alpha=0.7)
            gl.xlabels_top,gl.ylabels_right = False,False
            gl.xformatter,gl.yformatter = LONGITUDE_FORMATTER,LATITUDE_FORMATTER

            # add coastline, rivers, and lakes
            ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
            ax.add_feature(cfeature.RIVERS.with_scale('50m')) 
            ax.add_feature(cfeature.LAKES.with_scale('50m')) 
    
            # construct a meshgrid 
            XX, YY = np.meshgrid(lons_region,lats_region)
            Z = grid_region
        
            # calculate the maximum of the abs(Z)
            abs_Z_max = np.abs(Z).max()
            Z_levels = np.linspace(-abs_Z_max,abs_Z_max, 101)
        
            if block is None:
                CS = ax.contourf(XX,YY,Z,levels = Z_levels,extend='both',cmap=plt.cm.RdBu_r,transform = ccrs.PlateCarree())
            else:
                CS = ax.pcolormesh(gridlons,gridlats,Z,norm=colors.BoundaryNorm(boundaries=Z_levels,ncolors=256),cmap=plt.cm.RdBu_r,transform=ccrs.PlateCarree())

            #Create a colorbar and shrink it down a bit.
            cbar = plt.colorbar(CS,extend='both',format='%.0f',shrink=0.9)
            cbar.ax.set_ylabel(ylabel,fontsize=8)
            cbar.ax.tick_params(labelsize=8)
    
            if polygons is not None:
                ax.plot(polygons[:,1], polygons[:,0], color='m',transform=ccrs.Geodetic())

            if circles is not None:
                ax.scatter(circles[:,1], circles[:,0],facecolors="None", color='m',s=buffer_edge*30,transform=ccrs.Geodetic(),alpha=0.5)   
            
            plt.savefig(fig_dir+fig_name,bbox_inches='tight') 

    def set_region(self,region):
        info = self.info.copy()
        info['region'] = region
        lons_region,lats_region,grids_region,grids_std_region,lons_flag,lats_flag = crop_region(self.lons,self.lats,self.grids,self.grids_std,region)
        return Grid(info,grids_region,grids_std_region,lons_region,lats_region,lons_flag,lats_flag)

    def leakage_correction(self,method,r,nodes=None,study_area=None):

        info = self.info.copy()
        lats,lons,grids = self.lats,self.lons,self.grids
        lons_flag,lats_flag = self.lons_flag,self.lats_flag

        if method == 'forward_model':
            corrected_grids,corrected_grids_std = forward_model(grids,r)
            info['title'] = 'Signal leakage corrected(forward model method) ' + info['title']

        elif method == 'space_domain':  
            mas,mas_std = space_domain(self,nodes,study_area,r)
            corrected_grids = np.zeros_like(grids)
            corrected_grids_std = corrected_grids.copy()
            
            nodes_index = nodes.nodes_index
            for k in range(len(nodes_index)):
                i,j = nodes_index[k]
                corrected_grids[:,i,j] = mas[:,k]
                corrected_grids_std[:,i,j] = mas_std[:,k]
            info['title'] = 'Signal leakage corrected(space domain method) ' + info['title']   
               
        else:
            raise Exception('Only forward_model and space domain are avaliable for grids')   
        return Grid(info,corrected_grids,corrected_grids_std,lons,lats,lons_flag,lats_flag)    
         

    def __add__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_grids,self_grids_std = self.grids,self.grids_std
        other_grids,other_grids_std = other.grids,other.grids_std

        if 'rate' in self.title and 'rate' in other.title:
            add_grids = self_grids + other_grids
            add_grids_std = np.sqrt(self_grids_std**2 + other_grids_std**2)

        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('addition can not be completed between rate object and series object') 

        else:
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_grids = np.zeros((n,)+self.grids.shape[1:])
                assumed_other_grids = assumed_self_grids.copy()

                assumed_self_grids[self_existing_solution_flag] = self.grids
                assumed_other_grids[other_existing_solution_flag] = other.grids

                self_grids = assumed_self_grids[existing_solution_flag] 
                other_grids = assumed_other_grids[existing_solution_flag] 

                assumed_self_grids_std = np.zeros((n,)+self.grids_std.shape[1:])
                assumed_other_grids_std = assumed_self_grids_std.copy()

                assumed_self_grids_std[self_existing_solution_flag] = self.grids_std
                assumed_other_grids_std[other_existing_solution_flag] = other.grids_std

                self_grids_std = assumed_self_grids_std[existing_solution_flag] 
                other_grids_std = assumed_other_grids_std[existing_solution_flag] 

                solution_month = self.total_month[existing_solution_flag]
                solution_counts = len(solution_month)
                missing_solution_flag = ~existing_solution_flag
                missing_month = self.total_month[missing_solution_flag]
                missing_month_counts = len(missing_month)

                info['solution_month'] = solution_month
                info['solution_counts'] = solution_counts
                info['missing_month'] = missing_month
                info['missing_month_counts'] = missing_month_counts
                info['missing_solution_flag'] = missing_solution_flag

            add_grids = self_grids + other_grids
            add_grids_std = np.sqrt(self_grids_std**2 + other_grids_std**2)

        return Grid(info,add_grids,add_grids_std,self.lons,self.lats,self.lons_flag,self.lats_flag)
    
    def __sub__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_grids,self_grids_std = self.grids,self.grids_std
        other_grids,other_grids_std = other.grids,other.grids_std

        if 'rate' in self.title and 'rate' in other.title:
            sub_grids = self_grids - other_grids
            sub_grids_std = np.sqrt(self_grids_std**2 + other_grids_std**2)
        
        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('subtraction can not be completed between rate object and series object') 

        else:
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_grids = np.zeros((n,)+self.grids.shape[1:])
                assumed_other_grids = assumed_self_grids.copy()

                assumed_self_grids[self_existing_solution_flag] = self.grids
                assumed_other_grids[other_existing_solution_flag] = other.grids

                self_grids = assumed_self_grids[existing_solution_flag] 
                other_grids = assumed_other_grids[existing_solution_flag] 

                assumed_self_grids_std = np.zeros((n,)+self.grids_std.shape[1:])
                assumed_other_grids_std = assumed_self_grids_std.copy()

                assumed_self_grids_std[self_existing_solution_flag] = self.grids_std
                assumed_other_grids_std[other_existing_solution_flag] = other.grids_std

                self_grids_std = assumed_self_grids_std[existing_solution_flag] 
                other_grids_std = assumed_other_grids_std[existing_solution_flag] 

                solution_month = self.total_month[existing_solution_flag]
                solution_counts = len(solution_month)
                missing_solution_flag = ~existing_solution_flag
                missing_month = self.total_month[missing_solution_flag]
                missing_month_counts = len(missing_month)

                info['solution_month'] = solution_month
                info['solution_counts'] = solution_counts
                info['missing_month'] = missing_month
                info['missing_month_counts'] = missing_month_counts
                info['missing_solution_flag'] = missing_solution_flag
        
            sub_grids = self_grids - other_grids
            sub_grids_std = np.sqrt(self_grids_std**2 + other_grids_std**2)

        return Grid(info,sub_grids,sub_grids_std,self.lons,self.lats,self.lons_flag,self.lats_flag)
    
    def __mul__(self,other):
        pass
    
    def __truediv__(self,other):
        pass             
    