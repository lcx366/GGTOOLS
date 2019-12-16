import numpy as np

from .regular_gldas import regular_gldas

def lsm(lons,lats,data,q='twsc',mode='classic'):
    '''
    Calculate some land surface quantities from GLDAS grid data, such as Terrestrial Water Storage Changes(TWSC).
    This program provides two methods to estimate TWSC. The first one is the classic(traditional) technique, i.e. 'TWSC = SoilMoi + Accum_Snow + Canopy'.
    The second way is through the water balance equation, i.e. 'TWSC = precipitation - water_evaporation - surface_runoff - subsurface_runoff'. The difference in TWSC obtained by these two methods is minimal.
    Note: The TWSC here has been de-averaged. The unit of TWSC is [kg/m2] or [mm w.e.]

    Usage:
    twsc = lsm(lons,lats,data)
    twsc = lsm(lons,lats,data,mode='wbe')

    Inputs:
    lats -> [float array] latitudes of gldas grid data
    lons -> [float array] longitudes of gldas grid data
    data -> [float 2d array] gldas grids data

    Parameters:
    q -> [optional, str, default = 'twsc'] land surface quantities to be calculated from GLDAS grid data. Currently, only TWSC is avaliable.
    mode -> [optional, str, default = 'classic'] Method to eatimate TWSC. Avaliable options are 'classic' and 'wbe'. If 'classic', TWSC = SoilMoi + Accum_Snow + Canopy. If 'wbe', TWSC = precipitation - water_evaporation - surface_runoff - subsurface_runoff.
    
    Outputs:
    twsc -> [float 2d array] TWSC grids 
    '''
    info = {}
    info['normalization'] = 'fully normalized'
    info['permanent_tide'] = 'inclusive'
    info['earth_gravity_param'] = '3.9860044150E+14 m3/s2'
    info['mean_equator_radius'] = '6.3781363000E+06 m'
    info['background_gravity'] = 'GGM05C'
    info['unused_days'] = []
        
    month2seconds = 30*24*3600 # Number of seconds per month
    month2threehours = 30*8 # Number of 3-hour per month
    if q == 'twsc':
        twscs = []
        if mode == 'classic':
            for datum in data:
                # Get monthly soil moisture from surface to 2m deep, snow melt, and canopy water
                SoilMoi = datum['SoilMoi0_10cm_inst']+datum['SoilMoi10_40cm_inst']+datum['SoilMoi40_100cm_inst']+datum['SoilMoi100_200cm_inst']
                Accum_Snow = datum['SWE_inst']
                Canopy = datum['CanopInt_inst']
                # Calculate TWSC for each month
                twsc = SoilMoi + Accum_Snow + Canopy
                twscs.append(np.array(twsc[0]))
            twscs = np.array(twscs)
            # De-average TWSC
            twscs_deaverage = twscs - np.average(twscs,axis=0)
            info['summary'] = 'TWSC is estimated by the formulation: [TWSC = SoilMoi(0-200cm) + Accum_Snow + Canopy], and it has been converted to Equivalent Water Thickness in mm w.e.'
        
        elif mode == 'wbe': 
            for datum in data:         
                # Get monthly precipitation, evapotranspiration, surface runoff, and underground runoff
                precipitation_flux = datum['Rainf_f_tavg']
                water_evaporation_flux = datum['Evap_tavg']
                surface_runoff_amount = datum['Qs_acc']
                subsurface_runoff_amount = datum['Qsb_acc']
                # Calculate the TWSC for each month
                twsc = (precipitation_flux - water_evaporation_flux)*month2seconds - (surface_runoff_amount + subsurface_runoff_amount)*month2threehours
                twscs.append(np.array(twsc[0]))  
            twscs = np.array(twscs) 
            twscs_acc = np.cumsum(twscs,axis=0)

            # De-average TWSC
            twscs_deaverage = twscs_acc - np.average(twscs_acc,axis=0)
            info['summary'] = 'TWSC is estimated by the Water Balance Equation: [TWSC = precipitation - water_evaporation - surface_runoff - subsurface_runoff], and it has been converted to Equivalent Water Thickness in mm w.e.'

        else:
            raise Exception("Currenly, only 'classic' and 'wbe' are available for calculating TWSC.")  

        lats,lons,twscs_deaverage = regular_gldas(lats,lons,twscs_deaverage)  

        info['title'] =  'Terrestrial Water Storage Change(TWSC) derived from the '+ datum.title
        info['equi_material'] = 'Water' 
        info['processing_level'] = datum.conventions
        info['product_version'] = datum.source

        return info,lats,lons,twscs_deaverage

    else:
        raise Exception("Currenly, only 'twcs' is feasible.")