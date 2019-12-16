from os import path,makedirs
from urllib.request import urlretrieve

def static_download(model): 
    '''
    Download static gravity modle from icgem.gfz-potsdam.de; if the file to be downloaded is already included in the download directory, the download is automatically skipped.
    
    Usage: 
    static_download('GGM05C')
    static_download('EIGEN-6C4')
    
    Inputs:
    model -> [str] Available options are 'GGM05C' and 'EIGEN-6C'.
    
    Outputs: downloaded static gravity model
        
    Examples:
    >>> static_download('GGM05C')
    Downloading the static gravity model GGM05C ... Finished
    'static_models/GGM05C.gfc'
    >>> static_download('EIGEN-6C4')
    Downloading the static gravity model EIGEN-6C4 ... Finished
    'static_models/EIGEN-6C4.gfc'
    '''

    direc = 'static_models/'
    if not path.exists(direc): makedirs(direc)
    
    if model == 'GGM05C':
        gravity_file = direc + 'GGM05C.gfc'
        url = 'http://icgem.gfz-potsdam.de/getmodel/gfc/778a683780a5b0ad3163f4772b97b9075a0a13c389d2bd8ea3f891b64cfa383d/GGM05C.gfc'
    elif model == 'EIGEN-6C4':
        gravity_file = direc + 'EIGEN-6C4.gfc'
        url = 'http://icgem.gfz-potsdam.de/getmodel/gfc/7fd8fe44aa1518cd79ca84300aef4b41ddb2364aef9e82b7cdaabdb60a9053f1/EIGEN-6C4.gfc'
    else:
        raise Exception('Currently, available static gravity models are GGM05C and EIGEN-6C4.')    

    if not path.exists(gravity_file):
        print('Downloading the static gravity model '+ model,end=' ... ')
        urlretrieve(url, gravity_file)
        print('Finished')    
    return gravity_file    