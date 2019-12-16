from os import path,makedirs
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.signal import correlate

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import GPy

from ..gg.utils import month2int
from ..gg.fit_func import func
from ..gg.lsq import lsqm,ilsqm,wlsqm,iwlsqm

class Series(object):
    '''
    class Series
    - attributes:
        - info
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
        - total_month
        - total_month_counts
        - solution_month
        - solution_counts
        - missing_month
        - missing_month_counts
        - missing_solution_flag
        - unused_days
        - date_issued
        - equi_material 
        - filter
        - area
        - qs
        - qs_std
    
    - methods:
        - rate
        - plot
        - leakage_correction
    '''
    def __init__(self,info,area,qs,qs_std):
        
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
        self.solution_month = info['solution_month'] 
        self.solution_counts = info['solution_counts'] 
        self.missing_month = info['missing_month']
        self.missing_month_counts = info['missing_month_counts']
        self.missing_solution_flag = info['missing_solution_flag']
        self.unused_days = info['unused_days']
        self.date_issued = info['date_issued']
        self.equi_material = info['equi_material']
        self.filter = info['filter']
        
        self.area = area
        self.qs,self.qs_std = qs,qs_std
        
        
    def __repr__(self):
    
        return 'title = {:s}\ntime_coverage_start = {:s}\ntime_coverage_end = {:s}\nsolution_counts = {:d}\ntotal_month_counts = {:d}\nmissing_month_counts = {:d}\nequi_material = {:s}'.format\
        (self.title,self.time_coverage_start,self.time_coverage_end,self.solution_counts,self.total_month_counts,self.missing_month_counts,self.equi_material)
    
    def rate(self,mode='ILSQM'):
        '''
        Estimate the annual change rate of grids for Surface Mass Anomaly in Equivalent Water(or Ice, Sand) Thickness(EWT) using the linearly fitting method.
        There are four methods for linearly fitting. They are 'LSQM', 'ILSQM', 'WLSQM', and 'IWLSQM', respectively. 
        Note that 'LSQM' and 'WLSQM' will not remove any outliers caused by abnormal monthly solutions.
        Usage:
        xxx_sma_grid_rate = xxx_sma_grid.rate() or xxx_sma_grid_rate = xxx_sma_grid.rate('IWLSM')

        Parameters:
        lsm -> [optional, str, default = 'ILSQM'] Available options are 'LSM', 'ILSM', 'WLSM', and 'IWSLM', where
        'LSQM'   --  Least Square Method
        'ILSQM'  --  Iterative Least Square Method
        'WLSQM'  --  Weighted Least Square Method
        'IWLSQM' --  Iterative Weighted Least Square Method
        If the fitting method is not specified, it defaults to 'ILSQM'.
        If 'WLSQM' or 'IWLSQM' are used, the standard deviation for the time series will participate in the fitting process.
            
        Outputs:
        xxx_sma_grid_rate: instance of class Series
            
        Examples:
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
        qs,qs_std = self.qs,self.qs_std
        month = month2int(self.solution_month)
        qs_rate = qs_rate_std = 0
    
        if mode is 'LSQM':
            qs_rate,qs_rate_std,_,_, = lsqm(month,qs)
        elif mode is 'WLSQM':
            qs_rate,qs_rate_std,_,_ = wlsqm(month,qs,qs_std)
        elif mode is 'ILSQM':
            qs_rate,qs_rate_std,normal,_,_ = ilsqm(month,qs)
        elif mode is 'IWLSQM':
            qs_rate,qs_rate_std,normal,_,_ = iwlsqm(month,qs,qs_std)
        else:
            raise Exception('Currently, the least square method can only be LSQM, WLSQM, ILSQM, and IWLSQM.')
            
        info = self.info.copy()
        info['title'] = 'Annual change rate of ' + info['title']
        info['summary'] = info['summary'].replace('an estimate of','an estimate of annual change rate of')
        
        for em in ['w','i','s']:
            info['summary'] = info['summary'].replace('[mm '+em+'.e.]','[mm '+em+'.e./yr]')

        return Series(info,self.area,np.array([qs_rate]),np.array([qs_rate_std]))

    def plot(self,fig_name,ylabel,kernel='mat32+period',mode='ILSQM'):
        '''
        The time series consists of a linear trend, annual variation, and interannual variation. 
        The linear trend and interannual variation together make up the long-term trend. 
        The time series is fitted by the Gaussian Process Regression(GPR) technique. 
        Currently, available kernels for GPR here include 'rbf', 'rbf+period', 'mat32+period', 'mat52+period', 'linear+period', 'spline+period', and 'poly+period'. 
        If the kernel is not specified, it defaults to 'mat32+period'. 
        To extract the long-term trend from the GPR fitted curve, a Gaussian filter with a radius of 3 years is employed.
        '''

        if 'rate' in self.title:
            raise Exception('The shape of the series data to be plotted should contain multiple elements')

        fig_dir = 'figures/'
        if not path.exists(fig_dir): makedirs(fig_dir)     

        # 计算变化率
        qs,qs_std = self.qs,self.qs_std
        month = month2int(self.solution_month)
        qs_rate = qs_rate_std = 0
        normal = np.ones_like(qs,dtype=bool)
    
        if mode is 'LSQM':
            qs_rate,qs_rate_std,intercept,_, = lsqm(month,qs)
        elif mode is 'WLSQM':
            qs_rate,qs_rate_std,intercept,_ = wlsqm(month,qs,qs_std)
        elif mode is 'ILSQM':
            qs_rate,qs_rate_std,normal,intercept,_ = ilsqm(month,qs)
        elif mode is 'IWLSQM':
            qs_rate,qs_rate_std,normal,intercept,_ = iwlsqm(month,qs,qs_std)
        else:
            raise Exception('Currently, the least square method can only be LSQM, WLSQM, ILSQM, and IWLSQM.')
        
        total_month_index = np.arange(self.total_month_counts)
        fit_qs = func(total_month_index, intercept,qs_rate)    

        # GPR
        qs_deanomaly = qs[normal]    
        month_deanomaly = month[normal]
        # kernels
        if kernel == 'rbf':
            kernel = GPy.kern.RBF(input_dim=1)
        elif kernel == 'rbf+period':
            kernel = GPy.kern.RBF(input_dim=1)+GPy.kern.StdPeriodic(input_dim=1,period=12)
        elif kernel == 'mat32+period':
            kernel = GPy.kern.Matern32(input_dim=1)+GPy.kern.StdPeriodic(input_dim=1,period=12)
        elif kernel == 'mat52+period':
            kernel = GPy.kern.Matern32(input_dim=1)+GPy.kern.StdPeriodic(input_dim=1,period=12)    
        elif kernel == 'linear+period':
            kernel = GPy.kern.Linear(input_dim=1)+GPy.kern.StdPeriodic(input_dim=1,period=12)
        elif kernel == 'spline+period':
            kernel = GPy.kern.Spline(input_dim=1)+GPy.kern.StdPeriodic(input_dim=1,period=12)  
        elif kernel == 'poly+period':
            kernel = GPy.kern.Poly(input_dim=1)+GPy.kern.StdPeriodic(input_dim=1,period=12)             
        else:
            raise Exception('Currently, kernels can only be rbf, rbf+period, mat32+period, mat52+period, linear+period, spline+period, and poly+period.')
        
        model = GPy.models.GPRegression(month_deanomaly[:,None],qs_deanomaly[:,None],kernel)
        model.optimize()
        mean,var = model.predict(total_month_index[:,None])
        #print(model)

        month_date = np.array(self.total_month, dtype=np.datetime64)
        # plot

        plt.clf()
        fig = plt.figure(dpi=200)
        ax = fig.add_subplot(1, 1, 1)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

        plt.ylabel(ylabel)    
  
        plt.scatter(month_deanomaly,qs_deanomaly,marker = "+",s=25,label=None)
        plt.plot(total_month_index,mean[:,0],linewidth = 1.5,color='red',label=None)
        plt.plot(total_month_index,fit_qs,linewidth = 1.5,linestyle='dashed',color='black',label='Linear Trend')
        plt.plot(total_month_index,gaussian_filter(mean[:,0],9),linewidth = 1.5,color='green',label='Long-term Trend')
        plt.legend(loc='best', fontsize='x-small')

        plt.xticks(total_month_index[::24],month_date[::24])
        plt.setp(ax.get_xticklabels(),rotation=30)

        plt.savefig(fig_dir+fig_name,bbox_inches='tight')   

    def leakage_correction(self,method,k=1):
        '''
        The scale factor(or gain factor) is sensitive to the linear trend of the time series and not to the annual term. 
        Besides, it has local characteristics, i.e., different regions have different scale factors. 
        This method is only applicable when the SMA distribution derived from GRACE in the study area is similar to the TWSC distribution derived from GLDAS.
        '''

        info = self.info.copy()

        if method == 'scale_factor':
            info['title'] = 'Signal leakage corrected(Scale factor method) ' + info['title']
        else:
            raise Exception('Only scale_factor method is available')    

        return Series(info,self.area,k*self.qs,k*self.qs_std)         

    def __add__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_qs,self_qs_std = self.qs,self.qs_std
        other_qs,other_qs_std = other.qs,other.qs_std

        if 'rate' in self.title and 'rate' in other.title:
            add_qs = self_qs + other_qs
            add_qs_std = np.sqrt(self_qs_std**2 + other_qs_std**2)

        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('addition can not be completed between rate object and series object') 

        else:    
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_qs = np.zeros(n)
                assumed_other_qs = assumed_self_qs.copy()

                assumed_self_qs[self_existing_solution_flag] = self.qs
                assumed_other_qs[other_existing_solution_flag] = other.qs

                self_qs = assumed_self_qs[existing_solution_flag] 
                other_qs = assumed_other_qs[existing_solution_flag] 

                assumed_self_qs_std = np.zeros(n)
                assumed_other_qs_std = assumed_self_qs_std.copy()

                assumed_self_qs_std[self_existing_solution_flag] = self.qs_std
                assumed_other_qs_std[other_existing_solution_flag] = other.qs_std

                self_qs_std = assumed_self_qs_std[existing_solution_flag] 
                other_qs_std = assumed_other_qs_std[existing_solution_flag] 

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

            add_qs = self_qs + other_qs
            add_qs_std = np.sqrt(self_qs_std**2 + other_qs_std**2)

        return Series(info,self.area,add_qs,add_qs_std)
    
    def __sub__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_qs,self_qs_std = self.qs,self.qs_std
        other_qs,other_qs_std = other.qs,other.qs_std

        if 'rate' in self.title and 'rate' in other.title:
            sub_qs = self_qs - other_qs
            sub_qs_std = np.sqrt(self_qs_std**2 + other_qs_std**2)

        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('subtraction can not be completed between rate object and series object') 

        else:     
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_qs = np.zeros(n)
                assumed_other_qs = assumed_self_qs.copy()

                assumed_self_qs[self_existing_solution_flag] = self.qs
                assumed_other_qs[other_existing_solution_flag] = other.qs

                self_qs = assumed_self_qs[existing_solution_flag] 
                other_qs = assumed_other_qs[existing_solution_flag] 

                assumed_self_qs_std = np.zeros(n)
                assumed_other_qs_std = assumed_self_qs_std.copy()

                assumed_self_qs_std[self_existing_solution_flag] = self.qs_std
                assumed_other_qs_std[other_existing_solution_flag] = other.qs_std

                self_qs_std = assumed_self_qs_std[existing_solution_flag] 
                other_qs_std = assumed_other_qs_std[existing_solution_flag] 
           
            sub_qs = self_qs - other_qs
            sub_qs_std = np.sqrt(self_qs_std**2 + other_qs_std**2)

        return Series(info,self.area,sub_qs,sub_qs_std)
    
    def __mul__(self,other):
        pass
    
    def __truediv__(self,other):
        pass                 