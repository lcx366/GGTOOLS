import numpy as np
from os import path,remove,walk,makedirs
from ftplib import FTP
from astropy.time import Time

from ..ggclasses.class_SLR_C20 import SLR_C20

def slr_c20_download(RL = 'RL06'):
    '''
    Download SLR C20 data from isdcftp.gfz-potsdam.de; if the file to be downloaded is already included in the download directory, the download is automatically skipped.
    
    Usage: 
    SLR_C20_download()
    SLR_C20_download('RL05')
    SLR_C20_download('RL06')
    
    Parameters:
    RL -> [str] [optional, default = 'RL06'] release of the SLR C20 solutions; currently, only RL05 and RL06 are available
    
    Outputs: downloaded SLR C20 data
        
    Examples:
    >>> SLR_C20_download()
    Downloading ...  TN-11_C20_SLR_RL06.txt ... 226 Transfer complete
    '''
    
    if RL not in ['RL05','RL06']:
        raise Exception('Currently, only release RL05 and RL06 are feasible')
        
    server = 'isdcftp.gfz-potsdam.de'
    ftp = FTP(server) 
    ftp.login()
    dir_slr_c20_from = 'grace/DOCUMENTS/TECHNICAL_NOTES'

    ftp.cwd(dir_slr_c20_from)
    files_list = ftp.nlst()
    
    dir_slr_c20_to = 'SLR_' + RL
    
    if not path.exists(dir_slr_c20_to): makedirs(dir_slr_c20_to)   
    
    for file in files_list:
        if 'C20_SLR_'+RL in file: 
            slr_c20_file = path.join(dir_slr_c20_to,file)
            print('Downloading ... ',file,end=' ... ')
            
            while True:
                try:
                    with open(slr_c20_file, 'wb') as local_file:
                        res = ftp.retrbinary('RETR ' + file, local_file.write) # RETR is an FTP command
                        print(res)
                        local_file.close()   
                        break
                   
                except:
                    local_file.close()
                    remove(slr_c20_file)       
    ftp.quit()  
    ftp.close()
 
def read_slr_c20(start_date = None,end_date = None,RL = 'RL06'):
    
    '''
    Read the SLR C20 file. The SLR C20 file must be downloaded by the program SLR_C20_download() before calling this program.
    
    Usage: 
    slr_c20 = read_SLR_C20()
    slr_c20 = read_SLR_C20('2008-01')
    slr_c20 = read_SLR_C20('2008-01','2015-09')
    slr_c20 = read_SLR_C20(end_date = '2015-09')

    Parameters:
    start_date -> [optional, str, default = None] start date for the data to be read, for example, '2004-07'. If None, there is no limit on the start date.
    end_date -> [optional, str, default = None] end date for the data to be read, for example, '2004-07'. If None, there is no limit on the end date. 
    RL -> [optional, str, default = 'RL06'] release of the SLR C20 solutions; currently only RL06 is available

    Outputs:
    slr_c20 -> instance of SLR_C20 class
        
    Examples:
    >>> slr_c20 = read_SLR_C20()
    >>> print(slr_c20.background_gravity)
    GGM05C
    >>> print(slr_c20.mean_c20)
    -0.00048416945732
    >>> print(slr_c20.date_issued)
    Created September 4, 2019 - last month reported is August 2019.
    >>> print(slr_c20.title)
    GRACE Technical Note 11
    >>> print(slr_c20.institution)
    Center for Space Research, The University of Texas at Austin
    >>> print(slr_c20.summary)
    As a convenience to users who wish to use a replacement value for C20, a monthly C20 estimate time series is provided. These estimates are obtained from the analysis of Satellite Laser Ranging (SLR) data to five geodetic satellites: LAGEOS-1 and 2, Starlette, Stella and Ajisai. The background gravity satellites model used in the SLR analysis is consistent with the GRACE Release-06 processing, including the use of the same Atmosphere-Ocean De-aliasing product.
    >>> print(slr_c20.product_version)
    RL06
    >>> print(slr_c20.time_coverage_start)
    2002-04
    >>> print(slr_c20.time_coverage_end)
    2019-08
    >>> print(slr_c20.total_month_counts)
    209
    >>> print(slr_c20.solution_counts)
    176
    >>> print(slr_c20.missing_month)
    ['2002-06' '2002-07' '2003-06' '2011-01' '2011-06' '2012-05' '2012-10'
     '2013-03' '2013-08' '2013-09' '2014-02' '2014-07' '2014-12' '2015-06'
     '2015-10' '2015-11' '2016-04' '2016-09' '2016-10' '2017-02' '2017-07'
     '2017-08' '2017-09' '2017-10' '2017-11' '2017-12' '2018-01' '2018-02'
     '2018-03' '2018-04' '2018-05' '2018-08' '2018-09']
    >>> print(slr_c20.missing_month_counts)
    33
    >>> print(slr_c20.C20.shape,slr_c20.C20_DEMEAN.shape,slr_c20.C20_STD.shape)
    (176,) (176,) (176,)
    '''
    if RL != 'RL06':
        raise Exception('Currently, only release RL06 is feasible')
        
    slr_c20_dates,slr_c20_dates_interval = [],[]
    c20,c20_interval =[],[]
    c20_demean,c20_demean_interval = [],[]
    c20_std,c20_std_interval = [],[] 
    
    file_dir = 'SLR_'+RL
    for (dirname, dirs, files) in walk(file_dir): pass
    for file in files:
        if 'C20_SLR' in file:
            slr_c20_file = file
            
    info = {'normalization':'fully normalized','background_gravity':'GGM05C','product_version':RL}
    
    with open(path.join(dirname,slr_c20_file),'r',errors='ignore') as f:
        # read the header
        j = 0
        for line in f:
            if 'PRODUCT:' in line:
                break 
            
            words = [word.strip() for word in line.split(':') ]
                
            if words[0] == 'UPDATE HISTORY':
                info['date_issued'] = words[1]
                
            elif words[0] == 'NOTES':
                info['notes'] = words[1]
                
            elif words[0] == 'AUTHORS'   :
                info['institution'] = words[1]
                
            elif words[0] == 'DESCRIPTION':
                info['summary'] = words[1]     

            elif j == 2:
                info['title'] = words[0]
            elif j in range(7,9):
                info['notes'] = info['notes'] + ' ' + words[0] 
                
            elif j in range(33,35):
                info['institution'] = (info['institution'] + ', ' + words[0]).strip(', ')
                
            elif j in range(39,45):
                if j == 42: info['summary'] = info['summary'] + ' ' + line.strip()
                info['summary'] = (info['summary'] + ' ' + words[0]).strip()
                
            elif j == 65:
                info['mean_c20'] = float(words[1])
                
            else:
                pass
                
            j+=1  
            
        # read C20
        
        for line in f:
            words = line.replace('D', 'E').split()

            slr_c20_start_mjd, slr_c20_end_mjd = float(words[0]), float(words[-2])
            slr_c20_dates.append(Time(slr_c20_start_mjd,format = 'mjd').iso[:7])
            
            c20.append(words[2])
            c20_demean.append(words[3])
            c20_std.append(words[4])
    
    f.close()
    
    c20 = np.array(c20,dtype=np.float)
    c20_demean = np.array(c20_demean,dtype=np.float)*1e-10
    c20_std = np.array(c20_std,dtype=np.float)*1e-10
    
    num_solutions = len(slr_c20_dates)
        
    for i in range(num_solutions-1):
        if slr_c20_dates[i+1] == slr_c20_dates[i]:
            slr_c20_dates[i+1] = str(np.array(slr_c20_dates[i+1], dtype=np.datetime64)+ 1)  
                
    if start_date is None: start_date = slr_c20_dates[0]
    if end_date is None: end_date = slr_c20_dates[-1]  
            
    num_month = (int(end_date[:4]) - int(start_date[:4]))*12 + int(end_date[5:7]) - int(start_date[5:7]) + 1
    month_list = np.array(np.array(start_date, dtype=np.datetime64)+np.arange(num_month),dtype=np.str)
    
    for i in range(num_solutions):
        if slr_c20_dates[i] in month_list:
            slr_c20_dates_interval.append(slr_c20_dates[i])
            c20_interval.append(c20[i])
            c20_demean_interval.append(c20_demean[i])
            c20_std_interval.append(c20_std[i])
        
    solution_month = slr_c20_dates_interval
    solution_counts = len(solution_month)
    
    missing_solution_flag = ~np.in1d(month_list,solution_month)
    missing_month_list = month_list[missing_solution_flag]        
    missing_month_counts = len(missing_month_list) 
    
    c20,c20_demean,c20_std = np.array(c20_interval),np.array(c20_demean_interval),np.array(c20_std_interval)
    
    info['time_coverage_start'] = start_date
    info['time_coverage_end'] = end_date
    info['total_month'] = month_list
    info['total_month_counts'] = num_month
    info['solution_month'] = solution_month
    info['solution_counts'] = solution_counts
    info['missing_month'] = missing_month_list
    info['missing_month_counts'] = missing_month_counts
    info['missing_solution_flag'] = missing_solution_flag
        
    return SLR_C20(info,c20,c20_demean,c20_std)        