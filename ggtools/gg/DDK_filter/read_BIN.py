import sys
import struct
import numpy as np
from scipy.linalg import block_diag

def read_BIN(file,mode = 'packed'):
    '''
    Read the binary file containing symmetric/full or block diagonal matrices and associated vectors and parameters.

    Usage: 
    dat = read_BIN(file)
    dat = read_BIN(file, mode = 'packed')
    dat = read_BIN(file, mode = 'full')
    
    Inputs:
    file -> [str] Input file

    Parameters:
    mode -> [optional, str, default = 'packed'] Available options are 'packed' or 'full' form for the filter matrices.
    If 'packed', the matrix remains in packed form (dat['pack1'] field). If 'full', the matrix expands to its full form (dat['mat1'] field).
    Warning: the 'full' option may cause excessive RAM memory use with large matrices.

    Outputs:
    dat -> [dic]: Dictionary with the file content
    
    Notice: 
    This program is translated from the matlab/octave source code read_BIN.m written by Roelof Rietbroek 2016. 
    For more information, please refer to https://github.com/strawpants/GRACE-filter
    '''
    if mode == 'packed':
        unpack = False
    elif mode == 'full':
        unpack = True # unpack matrix in full size
    else:
        raise Exception("Only 'packed' or 'full' are avaliable.")
    
    # ckeck endian
    endian = sys.byteorder

    if endian == 'little':
        # open the binary file in little endian
        f = open(file,'rb')
    else:
        raise Exception('The endian of the binary file is little, but the endian of OS is big.')    
    
    dat = {}
    # read the data version and type from the binary file
    dat['version'] = f.read(8).decode().strip()
    dat['type'] = f.read(8).decode()
    dat['descr'] = f.read(80).decode().strip()

    for key in ['nints','ndbls','nval1','nval2']:
        dat[key] = struct.unpack('<I',f.read(4))[0]

    for key in ['pval1','pval2']:
        dat[key] = struct.unpack('<I',f.read(4))[0]

    dat['nvec'],dat['pval2'] = 0,1
    dat['nread'],dat['nval2'] = 0,dat['nval1']

    # read additional nblocks parameter
    nblocks = struct.unpack('<i',f.read(4))[0]

    lists = f.read(dat['nints']*24).decode().split()
    for element in lists:
        dat[element] = struct.unpack('<i',f.read(4))[0]
    
    lists = f.read(dat['ndbls']*24).decode().replace(':','').split()
    for element in lists:
        dat[element] = struct.unpack('<d',f.read(8))[0]
    
    # side description meta data
    lists = f.read(dat['nval1']*24).decode()
    dat['side1_d'] = [(lists[i:i+24]).replace('         ','') for i in range(0, len(lists), 24)] 

    # type specific meta data
    dat['blockind'] = np.array(struct.unpack('<'+str(nblocks)+'i',f.read(4*nblocks)))

    dat['side2_d'] = dat['side1_d']

    # read matrix data
    npack1 = dat['pval1']*dat['pval2']
    dat['pack1'] = np.array(struct.unpack('<'+str(npack1)+'d',f.read(8*npack1)))

    f.close() # close file

    if not unpack: return dat
    
    sz = dat['blockind'][0]
    dat['mat1'] = dat['pack1'][:sz**2].reshape(sz,sz).T

    shift1 = shift2 = sz**2

    for i in range(1,nblocks):
        sz = dat['blockind'][i] - dat['blockind'][i-1]
        shift2 = shift1 + sz**2
        dat['mat1'] = block_diag(dat['mat1'],dat['pack1'][shift1:shift2].reshape(sz,sz).T)
        shift1 = shift2
    del dat['pack1']

    return dat