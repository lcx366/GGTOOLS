import numpy as np

def filterSH(W,cilm,cilm_std = None):
    '''
    Filter spherical harmonic coefficients with a block diagonal fitler matrix.

    Usage:
    cilm_filter = filterSH(W,cilm)
    cilm_filter,cilm_std_filter = filterSH(W,cilm,cilm_std)

    Inputs:
    W -> [dic]: Dictionary containing the filter matrix (read by read_BIN.py)
    cilm -> [float 3d array] Spherical harmonic coefficients in matrix form. cilm = clm for i = 0; cilm = slm for i = 1

    Parameters:
    cilm_std -> [float 3d array] standard deviation for the spherical harmonic coefficients

    Outputs:
    cilm_filter -> [float 3d array] Filtered spherical harmonic coefficients
    cilm_std_filter -> [float 3d array] standard deviation for the filtered spherical harmonic coefficients

    Notice: 
    This program is translated from the matlab/octave source code filterSH.m written by Roelof Rietbroek 2016. 
    For more information, please refer to https://github.com/strawpants/GRACE-filter
    '''
    # Extract filter matrix
    # Maximum degree of the input coefficients
    lmax = cilm.shape[1] - 1

    # Extract the minimum and maximum degree supported by the filter matrix
    lmaxfilt,lminfilt = W['Lmax'],W['Lmin']    

    # Determine the output maximum degree (limited by either the filter or input data)
    lmaxout = min(lmax,lmaxfilt)

    # Reserve space for output (will have same size as input) and set to zero
    cilm_filter = np.zeros_like(cilm)   
    cilm_std_filter = np.zeros_like(cilm_std)

    # Loop parameter indicating the previous block number and the end position in the packed matrix of the previous block
    lastblckind,lastindex = 0,0

    # loop over the available blocks
    for iblk in range(W['Nblocks']):
        # Get the degree of the block from the block index
        degree = (iblk+1)//2
    
        # Break loop if the degrees of the block are larger than the degrees of the input
        if degree > lmaxout: break
        trig = (iblk + int(iblk > 0) + 1)%2      
    
        # Compute the size of the side of the stored block
        sz = W['blockind'][iblk] - lastblckind
    
        # Initialize the filter order block to a unit diagonal matrix
        blockn = np.identity(lmaxfilt+1-degree)
    
        # Minimum (stored) degree for this particular block (may be limited by the mininum degree supported by the filter)
        lminblk = max(lminfilt,degree)

        shift = lminblk - degree
        # unpack the stored filterblock (vector) in a fully occupied order block matrix
        blockn[shift:,shift:] = W['pack1'][lastindex:lastindex+sz**2].reshape(sz,sz).T
    
        # Filter the input coefficients (this is in fact just a matrix vector multiplication)
        if trig:
            cilm_filter[0,degree:lmaxout+1,degree] = np.dot(blockn[:lmaxout+1-degree,:lmaxout+1-degree],cilm[0,degree:lmaxout+1,degree])
        else:
            cilm_filter[1,degree:lmaxout+1,degree] = np.dot(blockn[:lmaxout+1-degree,:lmaxout+1-degree],cilm[1,degree:lmaxout+1,degree])
            
        if cilm_std is not None: 
            if trig:
                cilm_std_filter[0,degree:lmaxout+1,degree] = np.sqrt(np.dot(blockn[:lmaxout+1-degree,:lmaxout+1-degree]**2,cilm_std[0,degree:lmaxout+1,degree]**2))
            else:
                cilm_std_filter[1,degree:lmaxout+1,degree] = np.sqrt(np.dot(blockn[:lmaxout+1-degree,:lmaxout+1-degree]**2,cilm_std[1,degree:lmaxout+1,degree]**2))
            
        # Prepare the loop variables for next block
        lastblckind = W['blockind'][iblk]
        lastindex = lastindex + sz**2 

    if cilm_std is None: 
        return cilm_filter
    else:
        return cilm_filter,cilm_std_filter    

