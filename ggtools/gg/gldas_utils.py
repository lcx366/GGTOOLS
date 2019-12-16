import numpy as np
import xarray as xr
from scipy.interpolate import RectSphereBivariateSpline
from pyshtools.shclasses import SHCoeffs,SHGrid

from os import getenv,path,makedirs,walk,remove
import requests

from .utils import print_error_gldas
from ..ggclasses.class_GLDAS import GLDAS

def gldas_download(uid,passwd,start_date,end_date,res = '1deg',source = 'NOAH'):
    '''
    Download the GLDAS grid data over a period defined by the start date and end date and its documentation from urs.earthdata.nasa.gov
    If data to be downloaded are already in the download directory, the performance of download is automatically skipped.
    Currently, only data from NOAH are feasible. Avaliable resolutions are 1deg and 0.25deg.

    Usage:
    gldas_download(uid,passwd,start_date,end_date)
    gldas_download(uid,passwd,start_date,end_date,'0.25deg')

    Inputs:
    uid -> [str] username for logging into urs.earthdata.nasa.gov
    passwd -> [str] password for logging into urs.earthdata.nasa.gov
    start_date -> [str] start date for the data to be downloaded, for example, '2004-07'
    end_date -> [str] end date for the data to be downloaded, for example, '2013-11'
    
    Parameters:
    res -> [optional, default = '1deg'] Resolution of the GLDAS grid data. Avaliable options are '1deg' and '0.25deg'.
    source -> [optional, default = 'NOAH'] Publisher of the GLDAS grid data. Currently, only 'NOAH' are feasible.
 
    Outputs:
    GLDAS grid data stored in the GLDAS/NOAH10 or GLDAS/NOAH025 directory
        
    Examples:
    >>> uid,passwd = 'your_username','your_passwd'
    >>> start_date,end_date = '2002-07','2019-09'
    >>> gldas_download(uid,passwd,start_date,end_date)
    Downloading ...  GLDAS_NOAH10_M.A200207.021.nc4 ... 210 Transfer complete
    ...
    ...
    Downloading ...  GLDAS_NOAH10_M.A201909.021.nc4 ... 210 Transfer complete
    >>>
    >>> start_date,end_date = '2002-07','2019-09'
    >>> gldas_download(uid,passwd,start_date,end_date,'0.25deg')
    Downloading ...  GLDAS_NOAH025_M.A200207.021.nc4 ... 210 Transfer complete
    ...
    ...
    Downloading ...  GLDAS_NOAH025_M.A201909.021.nc4 ... 210 Transfer complete
    '''  
    val_res,dir_res = print_error_gldas(source, res)
        
    # Create a .netrc file in the home directory
    home = getenv('HOME')
    if not path.exists(home+'/.netrc'):
        netrc_file = open(home+'/.netrc','w')
        netrc_file.write('machine urs.earthdata.nasa.gov login '+uid+' password '+passwd)
        netrc_file.close()
    
    dir_gldas_to = 'GLDAS/'+source+dir_res
    if not path.exists(dir_gldas_to):
        makedirs(dir_gldas_to)

    num_month = (int(end_date[:4]) - int(start_date[:4]))*12 + int(end_date[5:7]) - int(start_date[5:7]) + 1
    month_list = list(np.array(np.array(start_date, dtype=np.datetime64)+np.arange(num_month),dtype=np.str))

    server = 'https://hydro1.gesdisc.eosdis.nasa.gov'
        
    for (dirname, dirs, files) in walk(dir_gldas_to): pass
    
    for month in month_list:
        ym = month.replace('-','')
        gldas_file = 'GLDAS_'+source+dir_res+'_M.A'+ym+'.021.nc4'
        if gldas_file not in files:
            print('Downloading ... ',gldas_file,end=' ... ')
            url = server+'/data/GLDAS/GLDAS_'+source+dir_res+'_M.2.1/'+ym[:4]+'/'+gldas_file
            result = requests.get(url)
            
            # If the download fails, try to download 3 times
            for idownload in range(3):
                try:
                    result.raise_for_status()
                    local_file = open(dir_gldas_to+'/'+gldas_file, 'wb')
                    local_file.write(result.content)
                    local_file.close()   
                    print(str(result.status_code)+' Transfer complete')
                    break
                except:
                    local_file.close()
                    remove(dir_gldas_to+'/'+gldas_file) 
            if idownload == 2: raise Exception('Server did not respond, file download failed')        
    
    # Download the documentation
    readme_file = 'README_GLDAS2.pdf'
    url = server+'/data/GLDAS/GLDAS_'+source+dir_res+'_M.2.1/doc/'+ readme_file                
    if not path.exists(dir_gldas_to+'/'+readme_file):
        result = requests.get(url)
        doc_file = open(dir_gldas_to+'/'+readme_file,'wb')
        doc_file.write(result.content)
        doc_file.close() 

def read_gldas(start_date=None,end_date=None,res = '1deg',source = 'NOAH'):
    '''
    Read the GLDAS files into a GLDAS class instance. Before calling this program, it is recommended to be able to download all the GLDAS data needed using the program gldas_download.
    
    Usage:
    gldas = read_gldas()
    gldas = read_gldas('2008-07','2018-12')
    gldas = read_gldas(start_date = '2010-01')
    gldas = read_gldas(end_date = '2018-12',res = '0.25deg')

    Inputs:
    
    Parameters:
    start_date: [optional, str, default = None] start date for the data to be read, for example, '2004-07'. If None, there is no limit on the start date.
    end_date: [optional, str, default = None] end date for the data to be read, for example, '2004-07'. If None, there is no limit on the end date. 
    If both start date and end_date are None, all files in the storage directory will be read.

    Outputs:
    gldas -> GLDAS class instance
        
    Examples: 
    >>> gldas = read_gldas()
    >>> print(gldas.title)
    GLDAS2.1 LIS land surface model output monthly mean
    >>> print(gldas.solution_counts)
    210
    >>> print(gldas.institution)
    NASA GSFC
    >>>
    >>> gldas = read_gldas('2008-07','2018-12')
    >>> print(gldas.time_coverage_start)
    2008-07
    >>> print(gldas.time_coverage_end)
    2018-12
    >>> print(gldas.solution_counts)
    126
    '''
    val_res,dir_res = print_error_gldas(source, res)
    
    # build empty lists to record the grid data infomation.
    filelist,filelist_interval = [],[]
    gldas_dates,gldas_dates_interval = [],[]
    date_issued = []
    data = []
    info = {}

    # go through all grid data files
    file_dir = 'GLDAS/'+source+dir_res
    for (dirname, dirs, files) in walk(file_dir): pass
    
    # sort these files by month sequences.
    files = np.sort(files)
    
    for filename in files:
        if filename.endswith('.nc4'):
            filelist.append(path.join(dirname,filename))
            gldas_dates.append(filename[-14:-10]+'-'+filename[-10:-8])

    num_solutions = len(gldas_dates) 
    
    if start_date is None: start_date = gldas_dates[0]
    if end_date is None: end_date = gldas_dates[-1]  

    num_month = (int(end_date[:4]) - int(start_date[:4]))*12 + int(end_date[5:7]) - int(start_date[5:7]) + 1
    month_list = np.array(np.array(start_date, dtype=np.datetime64)+np.arange(num_month),dtype=np.str)
    
    for i in range(num_solutions):
        if gldas_dates[i] in month_list:
            gldas_dates_interval.append(gldas_dates[i])
            filelist_interval.append(filelist[i])

    solution_month = gldas_dates_interval
    solution_counts = len(solution_month)
    missing_solution_flag = ~np.in1d(month_list,solution_month)
    missing_month_list = month_list[missing_solution_flag]
    missing_month_counts = len(missing_month_list)

    for filename in filelist_interval:
        datum = xr.open_dataset(filename)
        date_issued.append(datum.history[17:27])
        data.append(datum) 
    lons,lats = datum['lon'][:],datum['lat'][:]              
    
    info['title'] = datum.title
    info['summary'] = datum.comment
    info['resolution'] = res
    info['degree_order'] = int(90/float(res.strip('deg')))-1
    info['max_degree'] = info['max_order'] = info['degree_order']
    info['institution'] = datum.institution
    info['time_coverage_start'] = start_date
    info['time_coverage_end'] = end_date
    info['total_month'] = month_list
    info['total_month_counts'] = num_month
    info['solution_month'] = solution_month
    info['solution_counts'] = solution_counts
    info['missing_month'] = missing_month_list
    info['missing_month_counts'] = missing_month_counts
    info['missing_solution_flag'] = missing_solution_flag
    info['missing_value'] = datum.missing_value
    info['tavg'] = 'past 3-hour average'
    info['acc'] = 'past 3-hour accumulation'
    info['inst'] = 'instantaneous'
    info['date_issued'] = date_issued
    info['filter'] = 'none'

    return GLDAS(info,lons,lats,data)