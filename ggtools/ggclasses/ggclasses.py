from os import path,makedirs
import numpy as np
from scipy.ndimage import gaussian_filter

from pyshtools.shclasses import SHCoeffs,SHGrid,Slepian
from pyshtools.expand import spharm,SHExpandDH
from pyshtools.spectralanalysis import SlepianCoeffsToSH,Curve2Mask
from pyshtools.shio import read_icgem_gfc

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import GPy

from ..gg.static_models import static_download
from ..gg.lovenums import lovenums
from ..gg.utils import month2int,med,crop_region
from ..gg.fit_func import func
from ..gg.lsq import lsqm,ilsqm,wlsqm,iwlsqm
from ..gg.lsm import lsm
from ..gg.landmask import landmask
from ..gg.filter import filter_ddk,filter_gaussian,filter_gaussian_inverse
from ..gg.leakage import forward_model,space_domain,spectral_domain

class GSM(object):
    '''
    class GSM
    - attributes:
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
        - shc
        - shc_std
    
    - methods:
        - deaverage
        - debackground
        - replace_slr_c20
        - filter_ddk
        - filter_gaussian
        - sma
        - gsm
        - rate
        - grid
        - study_area
        - leakage_correction
    '''
    def __init__(self,info,shc,shc_std):
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
        
        self.shc = shc
        self.shc_std = shc_std
        
    def __repr__(self):
    
        return 'title = {:s}\nmax_degree = {:d}\nmax_order = {:d}\ndegree_order = {:d}\nnormalization = {:s}\ninstitution = {:s}\nprocessing_level = {:s}\nproduct_version = {:s}\ntime_coverage_start = {:s}\ntime_coverage_end = {:s}\nsolution_counts = {:d}\ntotal_month_counts = {:d}\nmissing_month_counts = {:d}'.format\
        (self.title,self.max_degree,self.max_order,self.degree_order,self.normalization,self.institution,self.processing_level,self.product_version,self.time_coverage_start,self.time_coverage_end,self.solution_counts,self.total_month_counts,self.missing_month_counts)
    
    def deaverage(self):
        
        '''
        Deaverage the GSM solutions of the GRACE and GRACE-FO RL06 products
    
        Usage: 
        xxx_gsm_d = xxx_gsm.deaverage()
   
        Outputs:
        xxx_gsm_d -> instance of SHC class
        
        Examples:
        >>> csr_gsm_d = csr_gsm.deaverage()
        >>> print(csr_gsm_d)
        ''' 
        info = self.info.copy()  

        shc_deaverage = self.shc - np.average(self.shc,axis=0)
        
        info['title'] = 'Deaveraged ' + info['title']
        info['background_gravity'] = 'Average of monthly solutions'
        return GSM(info,shc_deaverage,self.shc_std)

    def debackground(self):
        
        '''
        Debackground the GSM solutions of the GRACE and GRACE-FO RL06 products
    
        Usage: 
        xxx_gsm_d = xxx_gsm.debackground()
 
        Outputs:
        xxx_gsm_d -> instance of GSM class
        
        Examples:
        >>> csr_gsm_d = csr_gsm.debackground()
        >>> print(csr_gsm_d)
        ''' 
        info = self.info.copy()
        degree_order = self.degree_order
        background_gravity = self.background_gravity
        gravity_file = static_download(background_gravity)
        if background_gravity == 'GGM05C':
            cilm,gm,r0,errors = read_icgem_gfc(gravity_file ,lmax=degree_order,errors='calibrated')
        elif background_gravity == 'EIGEN-6C4':  
            cilm,gm,r0,errors = read_icgem_gfc(gravity_file,lmax=degree_order,errors='formal')
        else:
            raise Exception('Currently, available background gravity models are GGM05C and EIGEN-6C4')        

        shc_debackground = self.shc - cilm
        info['title'] = 'Debackgrounded ' + info['title']
        return GSM(info,shc_debackground,self.shc_std)        
    
    def replace_slr_c20(self,slr_c20):
        '''
        Replace the C20 values from the GSM files of the GRACE and GRACE-FO RL06 products with the 2nd degree terms from SLR measurements.
    
        Usage:
         xxx_gsm_r = xxx_gsm.replace_slr_c20(slr_c20)

        Inputs:
        slr_c20 -> instance of SLR_C20 class
            
        Outputs:
        xxx_gsm_r -> instance of GSM class
        
        Examples:
        >>> csr_gsm_r = csr_gsm.replace_slr_c20(slr_c20)
        >>> print(csr_gsm_r)
        '''
        shc,shc_std = self.shc.copy(),self.shc_std.copy()
        shc[:,0,2,0] = slr_c20.c20
        shc_std[:,0,2,0] = slr_c20.c20_std
        
        info = self.info.copy()
        info['title'] = info['title'] + ' with C20 replaced by the SLR measurements'
        info['summary'] = info['summary'] + ' Note that the 2nd-degree terms have been replaced with the values from SLR C20.' 
        
        return GSM(info,shc,shc_std)
    
    def filter_ddk(self,filter_type = 'DDK5'):
        
        '''
        Filt the deaveraged GSM SHC with the DDK filter, where DDK1,DDK2,...,DDK8 are avaliable.
    
        Usage:
        xxx_gsm_fddk = xxx_gsm.filt_DDK()

        Parameters:
        filt_type: [optional, str, default = 'DDK5'] types of DDK filter. Avaliable options are 'DDK1', 'DDK2',...,'DDK8'
            
        Outputs:
        xxx_gsm_fddk: instance of GSM class
        
        Examples:
        >>> slr_c20 = read_slr_c20(end_date='2017-06')
        >>> slr_c20_deaverage = slr_c20.deaverage()
        >>> csr_gsm = read_gsm('CSR',96,lmax=179,end_date='2017-06')
        >>> csr_gsm_d = csr_gsm.deaverage()
        >>> csr_gsm_r = csr_gsm_d.replace_slr_c20(slr_c20_deaverage)
        >>> csr_gsm_fddk = csr_gsm_r.filt_ddk('DDK5')
        >>> print(csr_gsm_fddk.title)
        >>> print(csr_gsm_fddk.summary)
        DDK5 filtered Deaveraged GRACE Geopotential Coefficients CSR RL06 with C20 replaced by the SLR measurements
        Spherical harmonic coefficients representing an estimate of the mean gravity field of Earth during the specified timespan derived from GRACE mission measurements. These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes. Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product. The 0th and 1st degree terms are excluded from CSR level-2. Note that the 2nd degree terms have been replaced by the C20 values from SLR. Also note that C20 values from SLR also experienced the DDK5 filtering.
        '''

        filter_shc,filter_shc_std = filter_ddk(filter_type,self.shc,self.shc_std)
            
        info = self.info.copy()
        info['title'] = filter_type + ' filtered ' + info['title']
        info['filter'] = filter_type
        
        if 'with C20 replaced by the SLR measurements' in info['title']:
            info['summary'] = info['summary'] + ' Also note that C20 from SLR also experienced the ' + filter_type + ' filtering.' 
          
        return GSM(info,filter_shc,filter_shc_std)
    
    def filter_gaussian(self,r):
        '''
        Filtering the deaveraged GSM SHC with the Gaussian filter.
    
        Usage: xxx_gsm_gau = xxx_gsm.filt_gaussian()

        Inputs:
        r -> [float] Gaussian filter radius in km
            
        Outputs:
        xxx_gsm_gau: instance of GSM class
        
        Examples:
        >>> slr_c20 = read_slr_c20(end_date='2017-06')
        >>> slr_c20_deaverage = slr_c20.deaverage()
        >>> csr_gsm = read_GSM('CSR',96,lmax=179,end_date='2017-06')
        >>> csr_gsm_d = csr_gsm.deaverage()
        >>> csr_gsm_r = csr_gsm_d.replace_slr_c20(slr_c20_deaverage)
        >>> csr_gsm_fgau = csr_gsm_r.filt_gaussian(200)
        >>> print(csr_gsm_fgau.title)
        Gaussian filtered Deaveraged GRACE Geopotential Coefficients CSR RL06 with C20 replaced by the SLR measurements
        >>> print(csr_gsm_fgau.summary)
        Spherical harmonic coefficients representing an estimate of the mean gravity field of Earth during the specified timespan derived from GRACE mission measurements. These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes. Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product. The 0th and 1st degree terms are excluded from CSR level-2. Note that the 2nd degree terms have been replaced by the C20 values from SLR. Also note that C20 values from SLR also experienced the Gaussian filtering.
        '''
        filter_shc,filter_shc_std = filter_gaussian(r,self.shc,self.shc_std)
        
        info = self.info.copy()
        info['title'] = 'Gaussian filtered ' + info['title']
        info['filter'] = 'Gaussian filter with radius of '+str(r) + ' km'
        
        if 'with C20 replaced by the SLR measurements' in info['title']:
            info['summary'] = info['summary'] + ' Also note that C20 from SLR also experienced the Gaussian filtering.' 
        
        return GSM(info,filter_shc,filter_shc_std)
    
    def sma(self, equi_material = None):
        '''
        Convert Stokes coefficents(or rates) for GSM to that for Surface Mass Anomaly in Equivalent Water(or Ice, Sand) Thickness(EWT) with unit of [mm w.e.](or [mm i.e.],[mm s.e.]) or [mm w.e./yr](or [mm i.e./yr],[mm s.e./yr]) 

        Usage: 
        xxx_sma = xxx_gsm.sma()

        Parameters:
        equi_material -> [optional, str, default = None] Equivalent material for Surface Mass Anomaly. Currently, only Water, Ice, and Sand are avaliable.

        Outputs:
        xxx_sma: instance of SMA class
            
        Examples:
        >>> csr_gsm = read_gsm('CSR',96)
        >>> gfz_gsm = read_gsm('GFZ',96)
        >>> jpl_gsm = read_gsm('JPL',96)
        >>> slr_c20 = read_slr_c20()
        >>> comb_gsm = GSM_average([csr_gsm.deaverage(),gfz_gsm.deaverage(),jpl_gsm.deaverage()])
        >>> comb_gsm_r = comb_gsm.replace_slr_c20(slr_c20.deaverage())
        >>> comb_gsm_ddk5 = comb_gsm_r.filt_DDK('DDK5')
        >>> sma = comb_gsm_ddk5.sma('Sand')
        >>> print(sma.title)
        Stokes coefficients for Surface Mass Anomaly(SMA) in Equivalent Sand Thickness(EWT) derived from the DDK5 filtered Combined Deaveraged GRACE & GRACE-FO Geopotential Coefficients CSR RL06, GFZ RL06, JPL RL06 with C20 replaced by the SLR measurements
        >>> print(sma.summary)
        Spherical harmonic coefficients representing an estimate of the Surface Mass Anomaly(SMA) expressed in terms of Equivalent Sand[1442kg/m3] Thickness(EWT) with unit of [mm s.e.] during the specified timespan derived from GRACE & GRACE-FO mission measurements. These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes. Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product. Note that the 2nd degree terms have been replaced by the C20 values from SLR. Also note that C20 values from SLR also experienced the DDK5 filtering.
        >>> print(sma.material)
        Sand
        '''
        if equi_material is None:
            equi_material = self.equi_material
            
        if equi_material == 'Water':
            rho = 1000
        elif equi_material == 'Ice':
            rho = 917
        elif equi_material == 'Sand':
            rho = 1442
        else:
            raise Exception('Currently, the equivalent material for SMA can only be Water, Ice, or Sand.')
            
        # Calculate the average density of the Earth
        G = 6.67430e-11
        GM = float(self.earth_gravity_param.partition('m3/s2')[0])
        a = float(self.mean_equator_radius.partition('m')[0])
        rho_ave = 3*GM/(4*G*np.pi*a**3)
        
        sma_shc = np.zeros_like(self.shc)
        sma_shc_std = np.zeros_like(self.shc_std)

        for l in range(self.degree_order+1):
            k_l =lovenums(l) 
            factor = a*rho_ave/(3*rho)*(2*l+1)/(1+k_l)
            #factor = a*h_l[l]/(1+k_l) # for vertical displacement
            sma_shc[:,:,l,:] = factor*self.shc[:,:,l,:]*1e3 # in mm
            sma_shc_std[:,:,l,:] = factor*self.shc_std[:,:,l,:]*1e3
            
        info = self.info.copy()
        if 'change rate' in info['title']:
            info['title'] = 'Stokes coefficients for annual change rate of Surface Mass Anomaly(SMA) in Equivalent ' + equi_material + ' Thickness(EWT) derived from the ' + info['title']
            info['summary'] = info['summary'].replace('mean gravity field of Earth','Surface Mass Anomaly(SMA) expressed in terms of Equivalent ' + equi_material + '['+ str(rho)+ 'kg/m3]' + ' Thickness(EWT) with unit of [mm '+equi_material[0].lower()+'.e./yr]')
        else:    
            info['title'] = 'Stokes coefficients for Surface Mass Anomaly(SMA) in Equivalent ' + equi_material + ' Thickness(EWT) derived from the ' + info['title']
            info['summary'] = info['summary'].replace('mean gravity field of Earth','Surface Mass Anomaly(SMA) expressed in terms of Equivalent ' + equi_material + '['+ str(rho)+ 'kg/m3]' + ' Thickness(EWT) with unit of [mm '+equi_material[0].lower()+'.e.]')
        info['equi_material'] = equi_material
        return GSM(info,sma_shc,sma_shc_std)
    
    def gsm(self):
        '''
        Convert Stokes coefficents for Surface Mass Anomaly in Equivalent Water(or Ice, Sand) Thickness(EWT) with unit of [mm w.e.] to that for GSM 

        Usage: xxx_gsm = xxx_sma.gsm()

        Parameters:
        -----------
            None
            
        Returns:
        -----------
            xxx_gsm: instance of GSM class
            
        Examples:
        -----------
        >>> csr_gsm = read_GSM('CSR',96)
        >>> gfz_gsm = read_GSM('GFZ',96)
        >>> jpl_gsm = read_GSM('JPL',96)
        >>> slr_c20 = read_SLR_C20()
        >>> comb_gsm = GSM_average([csr_gsm.deaverage(),gfz_gsm.deaverage(),jpl_gsm.deaverage()])
        >>> comb_gsm_r = comb_gsm.replace_slr_c20(slr_c20.deaverage())
        >>> comb_gsm_ddk5 = comb_gsm_r.filt_DDK('DDK5')
        >>> sma = comb_gsm_ddk5.sma('Sand')
        >>> gsm = sma.gsm()
        >>> print(sma.title)
        Stokes coefficients for Surface Mass Anomaly(SMA) in Equivalent Sand Thickness(EWT) derived from the DDK5 filtered Combined Deaveraged GRACE & GRACE-FO Geopotential Coefficients CSR RL06, GFZ RL06, JPL RL06 with C20 replaced by the SLR measurements
        >>> print(sma.summary)
        Spherical harmonic coefficients representing an estimate of the Surface Mass Anomaly(SMA) expressed in terms of Equivalent Sand[1442kg/m3] Thickness(EWT) with unit of [mm w.e.] during the specified timespan derived from GRACE & GRACE-FO mission measurements. These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes. Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product. Note that the 2nd degree terms have been replaced by the C20 values from SLR. Also note that C20 values from SLR also experienced the DDK5 filtering.
        >>> print(sma.material)
        Sand 
        >>> print(gsm.title)
        DDK5 filtered Combined Deaveraged GRACE & GRACE-FO Geopotential Coefficients CSR RL06, GFZ RL06, JPL RL06 with C20 replaced by the SLR measurements
        >>> print(gsm.summary)
        Spherical harmonic coefficients representing an estimate of the mean gravity field of Earth during the specified timespan derived from GRACE & GRACE-FO mission measurements. These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes. Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product. Note that the 2nd degree terms have been replaced by the C20 values from SLR. Also note that C20 values from SLR also experienced the DDK5 filtering.
        >>> print(comb_gsm_ddk5.SHC[100,0,30,20])
        1.8369657302246403e-12
        >>> print(sma.SHC[100,0,30,20])
        0.9456468755977168
        >>> print(gsm.SHC[100,0,30,20])
        1.8369657302246403e-12
        '''
        
        equi_material = self.equi_material
        
        if equi_material is 'Water':
            rho = 1000
        elif equi_material is 'Ice':
            rho = 917
        elif equi_material is 'Sand':
            rho = 1442
        else:
            raise Exception('Currently, the equivalent material can only be Water, Ice, or Sand.')
            
        # Calculate the average density of the Earth
        G = 6.67430e-11
        GM = float(self.earth_gravity_param.partition('m3/s2')[0])
        a = float(self.mean_equator_radius.partition('m')[0])
        rho_ave = 3*GM/(4*G*np.pi*a**3)
        
        gsm_shc = np.zeros_like(self.shc)
        gsm_shc_std = np.zeros_like(self.shc_std)
        
        for l in range(self.degree_order+1):
            k_l =lovenums(l) 
            factor = 3*rho/(a*rho_ave)*(1+k_l)/(2*l+1)
            gsm_shc[:,:,l,:] = factor*self.shc[:,:,l,:]/1e3 
            gsm_shc_std[:,:,l,:] = factor*self.shc_std[:,:,l,:]/1e3
            
        info = self.info.copy()
        if 'change rate' in info['title']:
            info['title'] = info['title'].replace('Stokes coefficients for annual change rate of Surface Mass Anomaly(SMA) in Equivalent ' + equi_material + ' Thickness(EWT) derived from the ','')
            info['title'] = info['title'].replace('Stokes coefficients for Surface Mass Anomaly(SMA) in Equivalent ' + equi_material + ' Thickness(EWT) derived from the ','')
            info['summary'] = info['summary'].replace('Surface Mass Anomaly(SMA) expressed in terms of Equivalent ' + equi_material + '['+ str(rho)+ 'kg/m3]' + ' Thickness(EWT) with unit of [mm '+equi_material[0].lower()+'.e./yr]','mean gravity field of Earth')
        else:    
            info['title'] = info['title'].replace('Stokes coefficients for Surface Mass Anomaly(SMA) in Equivalent ' + equi_material + ' Thickness(EWT) derived from the ','')
            info['summary'] = info['summary'].replace('Surface Mass Anomaly(SMA) expressed in terms of Equivalent ' + equi_material + '['+ str(rho)+ 'kg/m3]' + ' Thickness(EWT) with unit of [mm '+equi_material[0].lower()+'.e.]','mean gravity field of Earth')
        return GSM(info,gsm_shc,gsm_shc_std)
    
    def rate(self,mode='ILSQM'):
        '''
        Estimate the annual change rate of Geopotential coefficients or Stokes coefficents for Surface Mass Anomaly in Equivalent Water(or Ice, Sand) Thickness(EWT) using the linearly fitting method.
        There are four methods for linearly fitting, including 'LSM', 'ILSM', 'WLSM', and 'IWSLM'. The ILSM is default and recommended.
        
        Usage: xxx_sma_rate = xxx_sma.rate() or xxx_gsm_rate = xxx_gsm.rate('IWLSM')

        Parameters:
        -----------
            lsm [str] [optional, default: ILSM] alternertively, 'LSM', 'ILSM', 'WLSM', and 'IWSLM' are available, where
            'LSM'   --  Least Square Method
            'ILSM'  --  Iterative Least Square Method
            'WLSM'  --  Weighted Least Square Method
            'IWSLM' --  Iterative Weighted Least Square Method
            
        Returns:
        -----------
            xxx_gsm: instance of GSM class
            
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
        shc,shc_std = self.shc,self.shc_std
        month = month2int(self.solution_month)
        shc_rate,shc_rate_std = [np.zeros_like(shc[0]) for k in range(2)]
        degree = order = self.degree_order
        for i in range(2):
            for l in range(2,degree+1): # start with 1 if consider the motion of mass center 
                for m in range(order+1):
                    if i==1 and m==0 or m > l: continue
                    if mode is 'LSQM':
                        shc_rate[i,l,m],shc_rate_std[i,l,m],_,_ = lsqm(month,shc[:,i,l,m])
                    elif mode is 'WLSQM':
                        shc_rate[i,l,m],shc_rate_std[i,l,m],_,_ = wlsqm(month,shc[:,i,l,m],shc_std[:,i,l,m])
                    elif mode is 'ILSQM':
                        shc_rate[i,l,m],shc_rate_std[i,l,m],normal,_,_ = ilsqm(month,shc[:,i,l,m])
                    elif mode is 'IWLSQM':
                        shc_rate[i,l,m],shc_rate_std[i,l,m],normal,_,_ = iwlsqm(month,shc[:,i,l,m],shc_std[:,i,l,m])
                    else:
                        raise Exception('Currently, the least square method can only be LSQM, WLSQM, ILSQM, and IWLSQM.')
        info = self.info.copy()
        info['title'] = 'Annual change rate of ' + info['title']
        info['summary'] = info['summary'].replace('an estimate of','an estimate of annual change rate of')
        for em in ['w','i','s']:
            info['summary'] = info['summary'].replace('[mm '+em+'.e.]','[mm '+em+'.e./yr]')

        return GSM(info,np.array([shc_rate]),np.array([shc_rate_std]))
    
    def grid(self,region=None):
        '''
        Expand spherical harmonic coefficients into a regional grid
        Usage:
            xxx_sma_grid = xxx_sma.grid(region)
        Parameters:
            region: [float array] range of an area, for example, [96.0,120.0,21.0,39.0] means the boundary in lon and lat;
        
        Returns:
            an instance of GRID class
        '''
        # Create an empty list to contain the grid of EWT from GRACE for each month

        grids = []
        shcs,shcs_std = self.shc,self.shc_std  

        for cilm in shcs:
            coeffs_class = SHCoeffs.from_array(cilm)
            grids_class = coeffs_class.expand()
            grids.append(grids_class.data) 
        grids = np.array(grids)
        grids_std = np.zeros_like(grids)

        lons = grids_class.lons()
        lats = grids_class.lats()

        if region is not None:

            if 'rate' not in self.title:
                print('The calculation will take a few minutes, please be patient.')

            lons_region,lats_region,grids_region,grids_std_region,lons_flag,lats_flag = crop_region(lons,lats,grids,grids_std,region)
        
            # Convert SHCs_std to grids_std
            lmax = self.degree_order

            k = 0
            for shc_std in shcs_std:
                i = 0
                for theta in lats_region:
                    j = 0
                    for phi in lons_region:
                        ylm = spharm(lmax, 90-theta, phi)
                        grids_std_region[k,i,j] = np.sqrt(np.sum((ylm*shc_std)**2))
                        j+=1
                    i+=1 
                k+=1  
        else:
            region = 'globe'
            grids_region = grids
            grids_std_region = np.zeros_like(grids_region)
            lons_region,lats_region = lons,lats
            lons_flag = np.ones(len(lons_region),dtype=bool)
            lats_flag = np.ones(len(lats_region),dtype=bool)

            # Note: Since it takes a lot of time to calculate the uncertainties of the global grid data, the uncertainties are all set to zero.                   
        
        info = self.info.copy()
        info['title'] =  'Grids expanded from ' + info['title']
        info['summary'] = info['summary'].replace('Spherical harmonic coefficients','Grids')
        info['summary'] = info['summary'].replace('coefficients','grids')
        info['region'] = region
        
        return Grid(info,grids_region,grids_std_region,lons_region,lats_region,lons_flag,lats_flag)
    
    def study_area(self,points,north_pole=False,central_meridian=False):
        
        a = float(self.mean_equator_radius.partition('m')[0])/1e3 # km
        
        if self.equi_material is 'Water':
            rho = 1000
        elif self.equi_material is 'Ice':
            rho = 917
        elif self.equi_material is 'Sand':
            rho = 1442
       
        qs,qs_std = [],[]

        mask_grid = Curve2Mask(2*(self.degree_order+1),points,north_pole,sampling=2,centralmeridian=central_meridian)
        mask_shc = SHExpandDH(mask_grid,sampling=2)
        
        area = mask_shc[0,0,0]*4*np.pi*a**2 # km^2
        
        for shc in self.shc:
            q = np.sum(shc*mask_shc)*4*np.pi*a**2*rho/1e9 # Gt
            qs.append(q)
        for shc_std in self.shc_std:    
            q_std = np.sqrt(np.sum((shc_std*mask_shc)**2))*4*np.pi*a**2*rho/1e9 # Gt
            qs_std.append(q_std)
        qs,qs_std = np.array(qs),np.array(qs_std) 
        
        info = self.info.copy()
        info['title'] = 'Integral(over the study area) of '+info['title']
        
        return Series(info,area,qs,qs_std)

    def leakage_correction(self,method,r,nodes=None,nodes_index=None,study_area=None,mode=None,ratio=None):

        if method == 'filter_inverse':
            info = self.info.copy()
            corrected_shc,corrected_shc_std = filter_gaussian_inverse(r,self.shc,self.shc_std)
            info['title'] = 'Signal leakage corrected(filter inverse method) ' + info['title']
            return GSM(info,corrected_shc,corrected_shc_std)

        elif method == 'spectral_domain': 
            shcs = self.shc
            grids = self.grid() 
            info = grids.info
            lmax = self.degree_order
            n = (lmax+1)*2
            mas,mas_std = spectral_domain(shcs,nodes,study_area,r,mode,ratio)

            corrected_grids = np.zeros((shcs.shape[0],n,2*n))
            corrected_grids_std = corrected_grids.copy()

            for k in range(len(nodes_index)):
                i,j = nodes_index[k]
                corrected_grids[:,i,j] = mas[:,k]
                corrected_grids_std[:,i,j] = mas_std[:,k]
            info['title'] = 'Signal leakage corrected(spectral domain method) ' + info['title'] 

            return Grid(info,corrected_grids,corrected_grids_std,grids.lons,grids.lats,grids.lons_flag,grids.lats_flag)           
        else:
            raise Exception('Only filter inverse and spectral domain are avaliable')               

    def __add__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_shc,self_shc_std = self.shc,self.shc_std
        other_shc,other_shc_std = other.shc,other.shc_std

        if 'rate' in self.title and 'rate' in other.title:
            add_shc = self_shc + other_shc
            add_shc_std = np.sqrt(self_shc_std**2 + other_shc_std**2)

        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('Addition can not be completed between rate object and series object') 

        else:    
            if self.solution_counts != other.solution_counts: 
                raise Exception('The number of monthly solutions is not consistent between the two added objects.')
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag   

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_shc = np.zeros((n,)+self.shc.shape[1:])
                assumed_other_shc = assumed_self_shc.copy()

                assumed_self_shc[self_existing_solution_flag] = self.shc
                assumed_other_shc[other_existing_solution_flag] = other.shc

                self_shc = assumed_self_shc[existing_solution_flag] 
                other_shc = assumed_other_shc[existing_solution_flag] 

                assumed_self_shc_std = np.zeros((n,)+self.shc_std.shape[1:])
                assumed_other_shc_std = assumed_self_shc_std.copy()

                assumed_self_shc_std[self_existing_solution_flag] = self.shc_std
                assumed_other_shc_std[other_existing_solution_flag] = other.shc_std

                self_shc_std = assumed_self_shc_std[existing_solution_flag] 
                other_shc_std = assumed_other_shc_std[existing_solution_flag] 

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
        
            add_shc = self_shc + other_shc
            add_shc_std = np.sqrt(self_shc_std**2 + other_shc_std**2)
        
        return GSM(info,add_shc,add_shc_std)
    
    def __sub__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_shc,self_shc_std = self.shc,self.shc_std
        other_shc,other_shc_std = other.shc,other.shc_std

        if 'rate' in self.title and 'rate' in other.title:
            sub_shc = self_shc - other_shc
            sub_shc_std = np.sqrt(self_shc_std**2 + other_shc_std**2)

        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('Subtraction can not be completed between rate object and series object') 
        
        else:    
            if self.solution_counts != other.solution_counts: 
                raise Exception('The number of monthly solutions is not consistent between the two subtracted objects.')
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_shc = np.zeros((n,)+self.shc.shape[1:])
                assumed_other_shc = assumed_self_shc.copy()

                assumed_self_shc[self_existing_solution_flag] = self.shc
                assumed_other_shc[other_existing_solution_flag] = other.shc

                self_shc = assumed_self_shc[existing_solution_flag] 
                other_shc = assumed_other_shc[existing_solution_flag] 

                assumed_self_shc_std = np.zeros((n,)+self.shc_std.shape[1:])
                assumed_other_shc_std = assumed_self_shc_std.copy()

                assumed_self_shc_std[self_existing_solution_flag] = self.shc_std
                assumed_other_shc_std[other_existing_solution_flag] = other.shc_std

                self_shc_std = assumed_self_shc_std[existing_solution_flag] 
                other_shc_std = assumed_other_shc_std[existing_solution_flag] 

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
        
            sub_shc = self_shc - other_shc
            sub_shc_std = np.sqrt(self_shc_std**2 + other_shc_std**2)
        
        return GSM(info,sub_shc,sub_shc_std)
    
    def __mul__(self,other):
        pass
    
    def __truediv__(self,other):
        pass       

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

    def leakage_correction(self,method,r,nodes=None,nodes_index=None,study_area=None):

        info = self.info.copy()
        lats,lons,grids = self.lats,self.lons,self.grids
        lons_flag,lats_flag = self.lons_flag,self.lats_flag

        if method == 'forward_model':
            corrected_grids,corrected_grids_std = forward_model(grids,r)
            info['title'] = 'Signal leakage corrected(forward model method) ' + info['title']

        elif method == 'space_domain':  
            mas,mas_std = space_domain(lats_flag,lons_flag,grids,nodes,study_area,r)
            corrected_grids = np.zeros_like(grids)
            corrected_grids_std = corrected_grids.copy()

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
            if self.solution_counts != other.solution_counts: 
                raise Exception('The number of monthly solutions is not consistent between the two added objects.') 
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
            if self.solution_counts != other.solution_counts: 
                raise Exception('The number of monthly solutions is not consistent between the two subtracted objects.')
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
            if self.solution_counts != other.solution_counts: 
                raise Exception('The number of monthly solutions is not consistent between the two added objects.')
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
            if self.solution_counts != other.solution_counts: 
                raise Exception('The number of monthly solutions is not consistent between the two subtracted objects.')
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
        
#----------------------------------------------------------------------    
class SLR_C20(object):
    '''
    class SLR_C20
    - attributes:
        - info 
        - normalization
        - background_gravity
        - title
        - summary
        - institution
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
        - date_issued
        - mean_c20
        - c20
        - c20_demean
        - c20_std
    
    - methods:
        - deaverage
    '''
    def __init__(self,info,c20,c20_demean,c20_std):
        self.info = info
        self.normalization = info['normalization'] 
        self.background_gravity = info['background_gravity']
        self.title = info['title']
        self.summary = info['summary']
        self.institution = info['institution']
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
        self.date_issued = info['date_issued']
        self.mean_c20 = info['mean_c20']
        
        self.c20 = c20
        self.c20_demean = c20_demean
        self.c20_std = c20_std
        
    def __repr__(self):
    
        return 'title = {:s}\nnormalization = {:s}\ninstitution = {:s}\nproduct_version = {:s}\ntime_coverage_start = {:s}\ntime_coverage_end = {:s}\nsolution_counts = {:d}\ntotal_month_counts = {:d}\nmissing_month_counts = {:d}'.format\
        (self.title,self.normalization,self.institution,self.product_version,self.time_coverage_start,self.time_coverage_end,self.solution_counts,self.total_month_counts,self.missing_month_counts)
        
    def deaverage(self):
        
        '''
        Deaverage the SLR C20 data.
    
        Usage: 
        slr_c20_deaverage = slr_c20.deaverage()

        Outputs:
        slr_c20_deaverage -> Deaveraged SLR C20 solitions
        
        Examples:
        -----------
        >>> slr_c20 = read_SLR_C20_file()
        >>> slr_c20_deaverage = slr_c20.deaverage()
        >>> print(slr_c20_deaverage)
        [ 1.58283192e-10,  1.20785392e-10, -4.92720085e-11, -3.20084085e-11,
         -1.73326085e-11,  9.12899151e-12,  5.35746915e-11,  9.35680915e-11,
          8.24645915e-11,  1.23541892e-10,  1.03958091e-10,  9.27922915e-11,
          ...
          ...
         -1.44952909e-10, -1.42519708e-10, -9.05243085e-11, -7.39151085e-11,
         -6.21788085e-11, -2.70057085e-11, -4.02519085e-11,  1.70110915e-11,
          2.92669150e-12, -4.49227085e-11, -1.89951408e-10, -1.86842409e-10]
        '''
        info = self.info.copy()

        background_mean = np.average(self.c20)

        info['background_gravity'] = 'Average of monthly solutions'   
        info['mean_c20'] = background_mean
        c20_deaverage = self.c20 - background_mean
        c20_demean = c20_deaverage
        
        info['title'] = 'Deaveraged ' + info['title']
        
        return SLR_C20(info,c20_deaverage,c20_demean,self.c20_std)

    def __add__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_c20,self_c20_demean,self_c20_std = self.c20,self.c20_demean,self.c20_std 
        other_c20,other_c20_demean,other_c20_std = other.c20,other.c20_demean,other.c20_std

        if 'rate' in self.title and 'rate' in other.title:
            add_c20 = self_c20 + other_c20
            add_c20_demean = self_c20_demean + other_c20_demean
            add_c20_std = np.sqrt(self_c20_std**2 + other_c20_std**2)

        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('addition can not be completed between rate object and series object') 

        else: 
            if self.solution_counts != other.solution_counts: 
                raise Exception('The number of monthly solutions is not consistent between the two added objects.')
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_c20 = np.zeros(n)
                assumed_other_c20 = assumed_self_c20.copy()

                assumed_self_c20[self_existing_solution_flag] = self.c20
                assumed_other_c20[other_existing_solution_flag] = other.c20

                self_c20 = assumed_self_c20[existing_solution_flag] 
                other_c20 = assumed_other_c20[existing_solution_flag] 

                assumed_self_c20_demean = np.zeros(n)
                assumed_other_c20_demean = assumed_self_c20_demean.copy()

                assumed_self_c20_demean[self_existing_solution_flag] = self.c20_demean
                assumed_other_c20_demean[other_existing_solution_flag] = other.c20_demean

                self_c20_demean = assumed_self_c20_demean[existing_solution_flag] 
                other_c20_demean = assumed_other_C20_demean[existing_solution_flag] 

                assumed_self_c20_std = np.zeros(n)
                assumed_other_c20_std = assumed_self_c20_std.copy()

                assumed_self_c20_std[self_existing_solution_flag] = self.c20_std
                assumed_other_c20_std[other_existing_solution_flag] = other.c20_std

                self_c20_std = assumed_self_c20_std[existing_solution_flag] 
                other_c20_std = assumed_other_c20_std[existing_solution_flag] 

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
        
            add_c20 = self_c20 + other_c20
            add_c20_demean = self_c20_demean + other_c20_demean
            add_c20_std = np.sqrt(self_c20_std**2 + other_c20_std**2)

        return SLR_C20(info,add_c20,add_c20_demean,add_c20_std)
    
    def __sub__(self,other):

        info = self.info.copy()
        #info['title'] = 
        #info['summary'] =

        self_c20,self_c20_demean,self_c20_std = self.c20,self.c20_demean,self.c20_std
        other_c20,other_c20_demean,other_c20_std = other.c20,other.c20_demean,other.c20_std

        if 'rate' in self.title and 'rate' in other.title:
            sub_c20 = self_c20 - other_c20
            sub_c20_demean = self_c20_demean - other_c20_demean
            sub_c20_std = np.sqrt(self_c20_std**2 + other_c20_std**2)

        elif ('rate' in self.title and 'rate' not in other.title) or ('rate' not in self.title and 'rate' in other.title):
            raise Exception('subtraction can not be completed between rate object and series object') 

        else:     
            if self.solution_counts != other.solution_counts: 
                raise Exception('The number of monthly solutions is not consistent between the two subtracted objects.')
            self_existing_solution_flag = ~self.missing_solution_flag
            other_existing_solution_flag = ~other.missing_solution_flag

            if (self_existing_solution_flag != other_existing_solution_flag).any():
                existing_solution_flag = self_existing_solution_flag & other_existing_solution_flag

                n = len(existing_solution_flag)

                assumed_self_c20 = np.zeros(n)
                assumed_other_c20 = assumed_self_c20.copy()

                assumed_self_c20[self_existing_solution_flag] = self.c20
                assumed_other_c20[other_existing_solution_flag] = other.c20

                self_c20 = assumed_self_c20[existing_solution_flag] 
                other_c20 = assumed_other_c20[existing_solution_flag] 

                assumed_self_c20_demean = np.zeros(n)
                assumed_other_c20_demean = assumed_self_c20_demean.copy()

                assumed_self_c20_demean[self_existing_solution_flag] = self.c20_demean
                assumed_other_c20_demean[other_existing_solution_flag] = other.c20_demean

                self_c20_demean = assumed_self_c20_demean[existing_solution_flag] 
                other_c20_demean = assumed_other_c20_demean[existing_solution_flag] 

                assumed_self_c20_std = np.zeros(n)
                assumed_other_c20_std = assumed_self_d20_std.copy()

                assumed_self_c20_std[self_existing_solution_flag] = self.c20_std
                assumed_other_c20_std[other_existing_solution_flag] = other.c20_std

                self_c20_std = assumed_self_c20_std[existing_solution_flag] 
                other_c20_std = assumed_other_c20_std[existing_solution_flag] 

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

            sub_c20 = self_c20 - other_c20
            sub_c20_demean = self_c20_demean - other_c20_demean
            sub_c20_std = np.sqrt(self_c20_std**2 + other_c20_std**2)

        return SLR_C20(info,sub_c20,sub_c20_demean,sub_c20_std)
    
    def __mul__(self,other):
        pass
    
    def __truediv__(self,other):
        pass          

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
   
class Mascon(object):
    '''
    Mascon class
    '''
    def __init__(self,lons,lats,lons_span,lats_span,area_mascon,mascon_mmwe):
        self.lons = lons
        self.lats = lats
        self.lons_span = lons_span
        self.lats_span = lats_span
        self.area_mascon = area_mascon
        self.mascon_mmwe = mascon_mmwe
        
    def __repr__(self):
        return 'GRACE Mascon class instance'

    def rate(self,mode='ILSQM'):
        pass

    def study_area(self,points,north_pole=False,central_meridian=False):
        pass    

    def plot(self,fig_name,ylabel,polygons=None,circles=None):
        pass        

    def set_region(self,region):
        pass