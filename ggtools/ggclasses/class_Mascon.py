from os import path,makedirs
import numpy as np
   
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
        return 'instance of class Mascon'

    def rate(self,mode='ILSQM'):
        pass

    def study_area(self,points,north_pole=False,central_meridian=False):
        pass    

    def plot(self,fig_name,ylabel,polygons=None,circles=None):
        pass        

    def set_region(self,region):
        pass