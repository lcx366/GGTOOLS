import numpy as np
from os import path,remove,walk,makedirs
from ftplib import FTP
from datetime import date,timedelta
from gzip import GzipFile

from .utils import print_error_grace
from ..ggclasses.class_GSM import GSM

def parse_gsm_filename(gsm_filename):
    '''
    Parse GRACE GSM filenames.
    
    Usage: 
    gsm_date = parse_gsm_filename(GSM_filename)
    
    Inputs:
    gsm_filename -> [str] filename of the GRACE RL06 Level-2 solutions;

    Outputs:
    gsm_date -> [dictionary] start yy-mm-dd and end yy-mm-dd;
        
    Examples:
    >>> gsm_date = parse_gsm_filename('GSM-2_2009060-2009090_GRAC_UTCSR_BB01_0600')
    >>> print(gsm_date)
    {'start': {'year': '2009', 'month': '03', 'day': '01'}, 'end': {'year': '2009', 'month': '03', 'day': '31'}}
    '''
    # Convert [year, day of year] to [year, month, day]
    start_year,start_doy = int(gsm_filename[6:10]),int(gsm_filename[10:13])
    end_year,end_doy = int(gsm_filename[14:18]),int(gsm_filename[18:21])
    
    start_date = date(start_year, 1,1) + timedelta(days=start_doy-1)
    start_year = start_date.strftime('%Y')
    start_month = start_date.strftime('%m')
    start_day = start_date.strftime('%d')
    
    end_date = date(end_year, 1,1) + timedelta(days=end_doy-1)
    end_year = end_date.strftime('%Y')
    end_month = end_date.strftime('%m')
    end_day = end_date.strftime('%d')
    
    return {'start':{'year':start_year,'month':start_month,'day':start_day},\
            'end':{'year':end_year,'month':end_month,'day':end_day}} 
    
def print_gsm_date_coverage(source, D = 60, RL = 'RL06'):
    '''
    Print the date coverage for the GRACE GSM data from isdcftp.gfz-potsdam.de
    
    Usage: 
    print_gsm_date_coverage('CSR')
    print_gsm_date_coverage('GFZ',96)
    print_gsm_date_coverage('JPL',60)

    Inputs:
    source -> [str] source or publisher of the GRACE GSM Level-2 monthly solutions. Available options are 'CSR', 'GFZ', and 'JPL'.
    
    Parameters:
    D -> [optional, int, default = 60] degree of the GRACE Level-2 solutions; Available options are 60 or 96
    RL -> [optional, str, default = 'RL06'] release of the GRACE Level-2 solutions; currently only RL06 is available
    
    Outputs: messages for date coverage
        
    Examples:
    >>> print_gsm_date_coverage('CSR', 60, 'RL06')
    GSM data for CSR_RL06_60 is available from 2002-04-05 to 2017-06-28
    >>> print_gsm_date_coverage('GFZ')
    GSM data for GFZ_RL06_60 is available from 2002-04-05 to 2017-06-29
    >>> print_gsm_date_coverage('JPL', 96)
    GSM data for JPL_RL06_96 is available from 2002-04-04 to 2017-06-29
    '''
    print_error_grace(source, D, RL) 
    
    server = 'isdcftp.gfz-potsdam.de'
    ftp = FTP(server) 
    ftp.login() 
    dir_gsm_from1 = 'grace//Level-2/' + source + '/' + RL
    ftp.cwd(dir_gsm_from1)
    files_list1 = ftp.nlst()
    
    dir_gsm_from2 = '~/grace-fo/Level-2/' + source + '/' + RL
    ftp.cwd(dir_gsm_from2)
    files_list2 = ftp.nlst()
    
    files_list = files_list1 + files_list2
    
    ftp.quit()
    ftp.close()
    
    if D is 60:
        gsm_filenames = np.sort([s for s in files_list if 'GSM' in s and 'BA01' in s])
    if D is 96:
        gsm_filenames = np.sort([s for s in files_list if 'GSM' in s and 'BB01' in s])
        
    gsm_date = parse_gsm_filename(gsm_filenames[0])   
    start_year = gsm_date['start']['year']
    start_month = gsm_date['start']['month']
    start_day = gsm_date['start']['day']
    
    gsm_date = parse_gsm_filename(gsm_filenames[-1])
    end_year = gsm_date['end']['year']
    end_month = gsm_date['end']['month']
    end_day = gsm_date['end']['day']
    
    print('GSM data for {:s}_{:s}_{:d} is available from {:s}-{:s}-{:s} to {:s}-{:s}-{:s}'.format(source,RL,D,start_year,start_month,start_day,end_year,end_month,end_day))
    
def gsm_download(source,D = 60,start_date = None,end_date = None,RL = 'RL06'): 
    '''
    Download GRACE GSM data from isdcftp.gfz-potsdam.de; if the file to be downloaded is already included in the download directory, the download is automatically skipped.
    
    Usage: 
    gsm_download('CSR', 96)
    gsm_download('GFZ', 96,'2003-01','2016-07')
    gsm_download('JPL', 60,'2006-05')
    gsm_download('JPL', end_date = '2018-12')

    Inputs:
    source -> [str] Source or publisher of the GRACE GSM Level-2 monthly solutions. Available options are 'CSR', 'GFZ', and 'JPL'.
    
    Parameters:
    D -> [optional, int, default = 60] degree of the GRACE Level-2 solutions; Avaliable options are 60 or 96
    start_date -> [optional, str, default = None] start date of the data to be downloaded. If None, If None, there is no limit on the start date.
    end_date -> [optional, str, default = None] end date of the data to be downloaded. If None, there is no limit on the end date.
    RL -> [optional, str, default = 'RL06'] release of the GRACE Level-2 solutions; currently only RL06 is available
    
    Outputs: downloaded GRACE GSM data
        
    Examples:
    >>> gsm_download('CSR', 96)
    Downloading ...  GSM-2_2002095-2002120_GRAC_UTCSR_BB01_0600.gz ... 226 Transfer complete
    Downloading ...  GSM-2_2002123-2002137_GRAC_UTCSR_BB01_0600.gz ... 226 Transfer complete
    ...
    ...
    >>> gsm_download('JPL', 60,'2004-07','2013-06')
    >>> gsm_download('GFZ', 96,'2007-02')
    >>> gsm_download('GFZ', 60,end_date = '2015-11')
    '''
    print_error_grace(source, D, RL)
        
    server = 'isdcftp.gfz-potsdam.de'
    ftp = FTP(server) 
    ftp.login()
    dir_gsm_from1 = '~/grace//Level-2/' + source + '/' + RL + '/'
    ftp.cwd(dir_gsm_from1)
    files_list1 = ftp.nlst()
    
    dir_gsm_from2 = '~/grace-fo/Level-2/' + source + '/' + RL + '/'
    ftp.cwd(dir_gsm_from2)
    files_list2 = ftp.nlst()
    files_list = files_list1 + files_list2

    ftp.cwd('~/')
    
    if D is 60:
        gsm_filenames = np.sort([s for s in files_list if 'GSM' in s and 'BA01' in s])
    if D is 96:
        gsm_filenames = np.sort([s for s in files_list if 'GSM' in s and 'BB01' in s]) 
    
    if start_date is None:
        gsm_date = parse_gsm_filename(gsm_filenames[0])
        start_year = gsm_date['start']['year']
        start_month = gsm_date['start']['month']
        start_day = gsm_date['start']['day']
        start_date = date(int(start_year), int(start_month),int(start_day))
    else:
        start_date = date(int(start_date[:4]), int(start_date[5:7]),1)
        
    if end_date is None:
        gsm_date = parse_gsm_filename(gsm_filenames[-1])
        end_year = gsm_date['end']['year']
        end_month = gsm_date['end']['month']
        end_day = gsm_date['end']['day']
        end_date = date(int(end_year), int(end_month),int(end_day))

    else:
        if end_date[5:7] == '12':
            end_date = date(int(end_date[:4]), int(end_date[5:7]),31)
        else:
            end_date = date(int(end_date[:4]), int(end_date[5:7])+1,1) - timedelta(days=1)
    
    gsm_filenames_bydate = []
    for gsm_filename in gsm_filenames:
        gsm_date = parse_gsm_filename(gsm_filename) 
        
        gsm_start_year = gsm_date['start']['year']
        gsm_start_month = gsm_date['start']['month']
        gsm_start_day = gsm_date['start']['day']
        
        gsm_end_year = gsm_date['end']['year']
        gsm_end_month = gsm_date['end']['month']
        gsm_end_day = gsm_date['end']['day']
        
        gsm_start_date = date(int(gsm_start_year), int(gsm_start_month),int(gsm_start_day))
        gsm_end_date = date(int(gsm_end_year), int(gsm_end_month),int(gsm_end_day))
        
        if start_date < gsm_start_date < end_date or start_date < gsm_end_date < end_date:
            gsm_filenames_bydate.append(gsm_filename)
    
    dir_gsm_to = 'GRACE/'+source+'/'+RL+'_'+str(D)
    
    if not path.exists(dir_gsm_to): makedirs(dir_gsm_to)
        
    for (dirname, dirs, files) in walk(dir_gsm_to): pass
    
    for gsm_filename_bydate in gsm_filenames_bydate:
        if gsm_filename_bydate[:-3] not in files:
            gsm_file = path.join(dir_gsm_to,gsm_filename_bydate)
            
            print('Downloading ... ',gsm_filename_bydate,end=' ... ')
            
            # If the download fails, try to download 3 times
            for idownload in range(3):
                try:
                    with open(gsm_file, 'wb') as local_file:
                        if 'GRAC' in gsm_filename_bydate:
                            res = ftp.retrbinary('RETR ' + dir_gsm_from1 + gsm_filename_bydate, local_file.write) # RETR is an FTP command
                        if 'GRFO' in gsm_filename_bydate:
                            res = ftp.retrbinary('RETR ' + dir_gsm_from2 + gsm_filename_bydate, local_file.write) # RETR is an FTP command    
                        print(res)
                        local_file.close()   
                        break
                   
                except:
                    local_file.close()
                    remove(gsm_file)
            if idownload == 2: raise Exception('Server did not respond, file download failed')        

            g_file = GzipFile(gsm_file)
                
            open(gsm_file[:-3], "wb").write(g_file.read())
            g_file.close()
                
            remove(gsm_file)          
    
    ftp.quit()  
    ftp.close()

def parse_gsm_file(filename,lmax):
    '''
    Parse the GRACE GSM Level-2 file.

    Usage: 
    info,cilm,cilm_std = parse_gsm_file(filename,lmax)

    Inputs:
    filename -> [str] filename of the GRACE GSM solution;
    lmax -> [int] the maximum degree to read; 
        
    Outputs:
    info -> [dictionary] detailed information for GRACE GSM solutions
    cilm -> [float array] Dimensionless SHCs with shape of (2, lmax + 1, lmax + 1);
    cilm_std -> [float array] Dimensionless SHCs stds with shape of (2, lmax + 1, lmax + 1)
        
    Examples:
    >>> info,cilm,cilm_std = parse_gsm_file('GFZ/RL06_96/GSM-2_2016221-2016247_GRAC_GFZOP_BB01_0600',179)
    >>> info,cilm,cilm_std = parse_gsm_file('JPL/RL06_60/GSM-2_2003060-2003090_GRAC_JPLEM_BA01_0600',30) 
    '''
    info = {'degree_order':lmax}
    
    with open(filename,'r',errors='ignore') as f:
        # read the header
        for line in f:
            if '# End of YAML header' in line:
                break 
            
            words = [word.strip() for word in line.split(':') ]
        
            if words[0] == 'degree':
                info['max_degree'] = int(words[1])
                
            elif words[0] == 'order':
                info['max_order'] = int(words[1]) 
                
            elif words[0] == 'normalization':
                info['normalization'] = words[1]
              
            elif words[0] == 'permanent_tide_flag':
                info['permanent_tide'] = words[1].split()[0]    
                
            elif words[0] == 'value':
                if words[1] == '3.9860044150E+14' or words[1] == '3.9860044150e+14':
                    info['earth_gravity_param'] = words[1] + ' m3/s2'
                elif words[1] == '6.3781363000E+06' or words[1] == '6.3781363000e+06':
                    info['mean_equator_radius'] = words[1] + ' m'
                    info['background_gravity'] = 'GGM05C'
                elif words[1] == '6.3781364600E+06' or words[1] == '6.3781364600e+06':
                    info['mean_equator_radius'] = words[1] + ' m'
                    info['background_gravity'] = 'EIGEN-6C4'
                else:
                    raise Exception('unknown background gravity model')
                    
            elif words[0] == 'title':   
                info['title'] = words[1]
                
            elif words[0] == 'summary': 
                info['summary'] = words[1]
                
            elif words[0] == 'institution': 
                info['institution'] = words[1]  
                
            elif words[0] == 'processing_level': 
                info['processing_level'] = words[1]
                
            elif words[0] == 'product_version': 
                info['product_version'] = words[1]   
                
            elif words[0] == 'time_coverage_start': 
                info['time_coverage_start'] = words[1] + ':' + words[2] + ':' + words[3]      
            
            elif words[0] == 'time_coverage_end': 
                info['time_coverage_end'] = words[1] + ':' + words[2] + ':' + words[3]
                
            elif words[0] == 'unused_days':
                info['unused_days'] = words[1].strip('[ ]').split(', ')
                
            elif words[0] == 'date_issued': 
                info['date_issued'] = words[1] + ':' + words[2] + ':' + words[3] 
            else:
                pass

        # read SHCs
        
        cilm = np.zeros((2,lmax + 1, lmax + 1))
        cilm_std = np.zeros_like(cilm)
        cilm[0,0,0] = 1 # C00 = 1
        
        for line in f:
            words = line.replace('D', 'E').split()
            l, m = int(words[1]), int(words[2])
            
            if m > lmax: break
            if l > lmax: continue

            value_cs = [float(words[3]), float(words[4])]
            value_cs_std = [float(words[5]), float(words[6])]
            cilm[:, l, m] = value_cs
            cilm_std[:, l, m] = value_cs_std
    f.close()
    return info,cilm,cilm_std   

def read_gsm(source = 'CSR',D = 60,lmax = None,start_date = None,end_date = None,RL = 'RL06'):
    '''
    Read the GRACE GSM Level-2 files. Before calling this program, it is recommended to be able to download all the GRACE GSM data needed using the program GSM_download.
    
    Usage: 
    xxx_gsm = read_gsm(source,D,lmax,start_date,end_date)

    Parameters:
    source -> [optional, str, default = 'CSR'] source of the GRACE Level-2 solutions; vailable sources are CSR, GFZ, JPL
    D -> [optional, int, default = 60] degree of the GRACE Level-2 solutions. Avaliable options are 60 or 96
    lmax -> [optional, int, default = None] degree to be read; it should be any positive integer. if None, lmax = D
    start_date -> [optional, str, default = None] start date of the data to be read. If None, there is no limit on the start date.
    end_date -> [optional, str, default = None] end date of the data to be read. If None, there is no limit on the end date.
    RL -> [optional, str, default = 'RL06'] release of the GRACE Level-2 solutions; currently only RL06 is available

    Outputs:
    xxx_gsm -> instance of GSM class
        
    Examples:
    >>> csr_gsm = read_gsm('CSR',96)
    >>> print(CSR_GSM.degree_order)
    96
    >>> print(csr_gsm.max_degree,csr_gsm.max_order)
    96 96
    >>> print(csr_gsm.normalization)
    fully normalized
    >>> print(csr_gsm.earth_gravity_param)
    3.9860044150E+14 m3/s2
    >>> print(csr_gsm.mean_equator_radius)
    6.3781363000E+06 m
    >>> print(csr_gsm.background_gravity)
    GGM05C
    >>> print(csr_gsm.title)
    GRACE Geopotential Coefficients CSR RL06
    >>> print(csr_gsm.summary)
    Spherical harmonic coefficients representing an estimate of the mean gravity field of Earth during the specified timespan derived from GRACE mission measurements. These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes. Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product. The 0th and 1st degree terms are excluded from CSR level-2.
    >>> print(CSR_GSM.institution)
    UT-AUSTIN/CSR
    >>> print(csr_gsm.processing_level)
    2
    >>> print(csr_gsm.product_version)
    RL06
    >>> print(csr_gsm.time_coverage_start)
    2002-04
    >>> print(csr_gsm.time_coverage_end)
    2017-06
    >>> print(csr_gsm.solution_counts)
    163
    >>> print(csr_gsm.total_month_counts)
    183
    >>> print(csr_gsm.missing_month)
    ['2002-06' '2002-07' '2003-06' '2011-01' '2011-06' '2012-05' '2012-10'
     '2013-03' '2013-08' '2013-09' '2014-02' '2014-07' '2014-12' '2015-07'
     '2015-10' '2015-11' '2016-04' '2016-09' '2016-10' '2017-02']
    >>> print(CSR_GSM.missing_month_counts)
    20
    >>> print(csr_gsm.unused_days)
    ['2002-04-10', '2002-04-11', '2002-04-28', '2002-05-07', '2002-05-08', '2002-05-14', '2002-08-28', '2002-09-27', '2002-11-26', '2002-12-15', '2003-01-22', '2003-01-24', '2003-01-26', '2003-01-27', '2003-01-29', '2003-01-30', '2003-02-08', '2003-02-26', '2003-03-06', '2003-09-21', '2003-11-24', '2003-12-04', '2004-05-24', '2004-05-25', '2004-05-26', '2004-12-09', '2004-12-10', '2004-12-11', '2004-12-16', '2005-03-14', '2005-03-15', '2005-03-16', '2005-12-03', '2005-12-04', '2005-12-09', '2005-12-10', '2005-12-11', '2005-12-12', '2005-12-13', '2006-03-26', '2006-12-24', '2006-12-25', '2006-12-26', '2007-01-12', '2007-01-13', '2007-01-17', '2007-04-12', '2007-06-13', '2007-11-15', '2007-11-16', '2007-11-22', '2007-11-23', '2008-04-21', '2010-03-18', '2011-12-14', '2011-12-15', '2011-12-16', '2011-12-24', '2012-12-07', '2012-12-08', '2013-12-25', '2013-12-26', '2013-12-27', '2014-11-16', '2015-07-07', '2015-07-08', '2015-07-09', '2015-07-10', '2015-07-11', '2015-07-12', '2016-12-01', '2016-12-02', '2016-12-03', '2017-01-01', '2017-05-01', '2017-05-02', '2017-06-04']
    >>> print(CSR_GSM.SHC.shape,CSR_GSM.SHC_std.shape)
    >>> (163, 2, 97, 97) (163, 2, 97, 97)
    >>> gfz_gsm = read_gsm('GFZ',96, 179, '2007-05','2012-05')
    >>> print(gfz_gsm.info)
    {'degree_order': 179, 'max_degree': 96, 'max_order': 96, 'normalization': 'fully normalized', 'permanent_tide': 'exclusive', 'earth_gravity_param': '3.9860044150E+14 m3/s2', 'mean_equator_radius': '6.3781364600E+06 m', 'background_gravity': 'EIGEN-6C4', 'title': 'GRACE Geopotential GSM Coefficients GFZ RL06', 'summary': "Spherical harmonic coefficients representing an estimate of Earth's mean gravity field during the specified timespan derived from GRACE mission measurements. These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes. Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product.", 'institution': 'GFZ German Research Centre for Geosciences', 'processing_level': '2', 'product_version': '6.0', 'time_coverage_start': '2007-05', 'time_coverage_end': '2012-05', 'unused_days': ['2007-06-13', '2007-11-15', '2007-11-22', '2007-11-23', '2008-04-21', '2008-04-22', '2010-03-18', '2011-12-14', '2011-12-15', '2011-12-16'], 'date_issued': ['2018-09-26', '2018-09-27', '2018-08-02'], 'solution_counts': 58, 'total_month_counts': 61, 'missing_month': array(['2011-01', '2011-06', '2012-05'], dtype='<U7'), 'missing_month_counts': 3}
    >>> print(gfz_gsm.SHC.shape,gfz_gsm.SHC_std.shape)
    (58, 2, 180, 180) (58, 2, 180, 180)
    >>> jpl_gsm = read_gsm('JPL',60,30,'2010-01')
    >>> print(jpl_gsm.info)
    {'degree_order': 30, 'max_degree': 60, 'max_order': 60, 'normalization': 'fully normalized', 'permanent_tide': 'inclusive', 'earth_gravity_param': '3.9860044150e+14 m3/s2', 'mean_equator_radius': '6.3781363000e+06 m', 'background_gravity': 'GGM05C', 'title': 'GRACE Geopotential Coefficients JPL RL06', 'summary': "Spherical harmonic coefficients representing an estimate of Earth's mean gravity field during the specified timespan derived from GRACE mission measurements.  These coefficients represent the full magnitude of land hydrology, ice, and solid Earth processes.  Further, they represent atmospheric and oceanic processes not captured in the accompanying GAC product.", 'institution': 'NASA/JPL', 'processing_level': '2', 'product_version': '6.0', 'time_coverage_start': '2010-01', 'time_coverage_end': '2017-06', 'unused_days': ['2010-06-15', '2010-06-17', '2010-06-18', '2010-06-19', '2010-06-20', '2010-06-21', '2010-06-22', '2010-06-23', '2012-12-07', '2012-12-08', '2013-06-12', '2013-12-25', '2013-12-26', '2013-12-27', '2015-07-07', '2015-07-08', '2015-07-09', '2015-07-10', '2015-07-11', '2015-07-12', '2015-12-31', '2016-01-20', '2016-05-25'], 'date_issued': ['2018-05-20', '2018-08-19'], 'solution_counts': 73, 'total_month_counts': 90, 'missing_month': array(['2011-01', '2011-06', '2012-05', '2012-10', '2013-03', '2013-08',
     '2013-09', '2014-02', '2014-07', '2014-12', '2015-07', '2015-10',
     '2015-11', '2016-04', '2016-09', '2016-10', '2017-02'], dtype='<U7'), 'missing_month_counts': 17}
    >>> print(jpl_gsm.SHC.shape,jpl_gsm.SHC_std.shape)
    (73, 2, 31, 31) (73, 2, 31, 31)  
    ''' 
    print_error_grace(source, D, RL)
    
    if lmax is None: lmax = D

    # record the SHCs for each month.
    filelist,filelist_interval = [],[]
    gsm_dates,gsm_dates_interval = [],[]
    shc,shc_std = [],[]
    unused_days,date_issued = [],[]
    
    # Go through all GRACE GSM files and put them into the file list.
    file_dir = 'GRACE/'+source+'/'+RL+'_'+str(D)
    for (dirname, dirs, files) in walk(file_dir): pass

    # Sort files by month sequences.
    files = np.sort(files)
    
    for filename in files:
        if 'GSM' in filename:
            filelist.append(path.join(dirname,filename)) 
            gsm_date = parse_gsm_filename(filename) 
            gsm_dates.append(gsm_date['start']['year']+'-'+gsm_date['start']['month'])
                
    num_solutions = len(gsm_dates)
    
    for i in range(num_solutions-1):
        if gsm_dates[i+1] == gsm_dates[i]:
            gsm_dates[i+1] = str(np.array(gsm_dates[i+1], dtype=np.datetime64)+ 1)      
    
    if start_date is None: start_date = gsm_dates[0]
    if end_date is None: end_date = gsm_dates[-1]  
            
    num_month = (int(end_date[:4]) - int(start_date[:4]))*12 + int(end_date[5:7]) - int(start_date[5:7]) + 1
    month_list = np.array(np.array(start_date, dtype=np.datetime64)+np.arange(num_month),dtype=np.str)
    
    for i in range(num_solutions):
        if gsm_dates[i] in month_list:
            gsm_dates_interval.append(gsm_dates[i])
            filelist_interval.append(filelist[i])

    solution_month = gsm_dates_interval
    solution_counts = len(solution_month)
    missing_solution_flag = ~np.in1d(month_list,solution_month)
    missing_month_list = month_list[missing_solution_flag]        
    missing_month_counts = len(missing_month_list)  
    
    for filename in filelist_interval:
        info,cilm,cilm_std = parse_gsm_file(filename,lmax)
        #unused_days += info['unused_days']
        unused_days.append(info['unused_days'])
        date_issued.append(info['date_issued'][:10])
        shc.append(cilm)
        shc_std.append(cilm_std)

    shc,shc_std = np.array(shc),np.array(shc_std)         
    
    info['time_coverage_start'] = start_date
    info['time_coverage_end'] = end_date

    info['total_month'] = month_list
    info['total_month_counts'] = num_month

    info['solution_month'] = solution_month
    info['solution_counts'] = solution_counts
    
    
    info['missing_month'] = missing_month_list
    info['missing_month_counts'] = missing_month_counts

    info['missing_solution_flag'] = missing_solution_flag

    #info['unused_days'] = list(filter(''.__ne__, unused_days))
    info['unused_days'] = unused_days
    info['date_issued'] = date_issued
    info['title'] = info['title'].replace('GRACE-FO','GRACE & GRACE-FO')
    info['summary'] = info['summary'].replace('GRACE-FO','GRACE & GRACE-FO')
    info['equi_material'] = 'Water'
    info['filter'] = 'none'
           
    return GSM(info,shc,shc_std)

def gsm_average(gsm_list):
    '''
    Combine the (deaveraged) GSM solution from multiple institutions into an integrated one. The combined solution
    is defined as the average of these solutions.
    
    Usage: 
    comb_gsm = GSM_average([csr_gsm,gfz_gsm,jpl_gsm])

    Inputs:
    gsm_list -> [str list] The list of instance of GSM class for multiple institutions.
            
    Outputs:
    comb_gsm -> the instance of GSM class for combined GSM solutions
        
    Examples:
    >>> csr_gsm = read_gsm('CSR',96,start_date='2002-05',end_date='2019-07')
    >>> gfz_gsm = read_gsm('GFZ',96,end_date='2019-06')
    >>> jpl_gsm = read_gsm('JPL',96,start_date='2002-06')
    >>> comb_gsm = gsm_average([csr_gsm.deaverage(),gfz_gsm.deaverage(),jpl_gsm.deaverage()])
    >>> print(comb_gsm.title)
    Combined Deaveraged GRACE & GRACE-FO Geopotential Coefficients CSR RL06, GFZ RL06, JPL RL06
    >>> print(comb_gsm.institution)
    UT-AUSTIN/CSR, GFZ German Research Centre for Geosciences, NASA/JPL
    '''
    n = len(gsm_list)
    dic_shc,dic_shc_std = [],[]
    shc,shc_std = [],[]
    gsm_keys = {}.keys()
    
    for i in range(1,n):
        if gsm_list[i].degree_order != gsm_list[0].degree_order:
            raise Exception('Degree and order for gsm in gsm list are not identical.')
    dim = gsm_list[0].degree_order+1
    
    for j in range(n):
        # convert Stokes coefficients for a gravity model that is not GGM05C to those for GGM05C
        if gsm_list[j].mean_equator_radius == '6.3781364600E+06 m':
            ratio_r = float(gsm_list[j].mean_equator_radius.partition('m')[0])/6.3781363000E6
            ratio_rl = np.array([ratio_r**l for l in range(dim)])[:,None]
            ratio_gm = float(gsm_list[j].earth_gravity_param.partition('m3/s2')[0])/3.9860044150E14
            dic_shc.append(dict(zip(gsm_list[j].solution_month,gsm_list[j].shc*ratio_rl*ratio_gm))) 
            dic_shc_std.append(dict(zip(gsm_list[j].solution_month,gsm_list[j].shc_std*ratio_rl*ratio_gm)))
                      
        dic_shc.append(dict(zip(gsm_list[j].solution_month,gsm_list[j].shc)))
        dic_shc_std.append(dict(zip(gsm_list[j].solution_month,gsm_list[j].shc_std)))
        gsm_keys = gsm_keys | dic_shc[j].keys()
    gsm_keys = np.sort(list(gsm_keys))
   
    nan_matrix = np.full((2,dim,dim),np.nan)
    none_matrix = np.full((2,dim,dim),None)

    for key in gsm_keys:
        temp1,temp2,temp3 = [],[],[]
        
        for k in range(n):
            temp1.append(dic_shc[k].get(key,nan_matrix))
            temp2.append((dic_shc_std[k]).get(key,nan_matrix)**2)
            temp3.append((dic_shc[k].get(key) != none_matrix).any())
        m = np.count_nonzero(temp3)    

        shc.append(np.nanmean(temp1,axis=0))
        shc_std.append(np.sqrt(np.nanmean(temp2,axis=0)/m))
    shc,shc_std = np.array(shc),np.array(shc_std)   
    
    start_date = np.sort([gsm_list[k].info['time_coverage_start'] for k in range(n)])[0]
    end_date = np.sort([gsm_list[k].info['time_coverage_end'] for k in range(n)])[-1]
    num_month = (int(end_date[:4]) - int(start_date[:4]))*12 + int(end_date[5:7]) - int(start_date[5:7]) + 1
    month_list = np.array(np.array(start_date, dtype=np.datetime64)+np.arange(num_month),dtype=np.str)
    solution_month = list(gsm_keys)
    solution_counts = len(solution_month)
    missing_solution_flag = ~np.in1d(month_list,solution_month)
    missing_month_list = month_list[missing_solution_flag]       
    missing_month_counts = len(missing_month_list) 
    
    info = gsm_list[0].info.copy()
    title = gsm_list[0].info['title']
    info['permanent_tide'] = 'inclusive'

    if gsm_list[0].background_gravity == 'Average of monthly solutions':
        info['background_gravity'] = 'Average of monthly solutions'
    else:
        info['background_gravity'] = 'GGM05C'
    info['earth_gravity_param'] = '3.9860044150E+14 m3/s2'
    info['mean_equator_radius'] = '6.3781363000E+06 m'

    info['title'] = 'Combined ' + title.replace(title[-9:-5],'')
    info['summary'] = info['summary'].partition("product.")[0] + 'product.'
    info['time_coverage_start'] = start_date
    info['time_coverage_end'] = end_date

    info['total_month'] = month_list
    info['total_month_counts'] = num_month

    info['solution_month'] = solution_month
    info['solution_counts'] = solution_counts

    info['missing_month'] = missing_month_list
    info['missing_month_counts'] = missing_month_counts

    info['missing_solution_flag'] = missing_solution_flag
    
    info['unused_days'] = 'Invalid'
    info['date_issued'] = 'Invalid'

    for k in range(1,n):
        info['institution'] = info['institution'] + ', ' + gsm_list[k].info['institution']
    return GSM(info,shc,shc_std)