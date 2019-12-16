import numpy as np
from scipy.interpolate import interp1d

def lovenums(l):
    '''
    Estimate Load Love Numbers(LLNs). Intermediate numbers can be linearly interpolated with errors of less than 
    0.05% for all l<200, as opposed to direct calculation. Note: the first-degree LLNs is taken as 0.021.

    Usage: 
    k_l =lovenums(l)

    Inputs:
    l -> [int, int list or array] Degree of LLNs 

    Outputs:
    k_l -> [float, float array] LLNs for degree l
        
    Examples:
    >>> k_l = lovenums(45)
    >>> print(k_l)
    -0.03
    >>> k_ls = lovenums([45,58,96])
    >>> print(k_ls)
    [-0.03   -0.0242 -0.0148]
    >>> k_ls = lovenums(np.arange(10,20))
    >>> print(k_ls)
    [-0.069  -0.0665 -0.064  -0.062  -0.06   -0.058  -0.0566 -0.0552 -0.0538 -0.0524]
    ''' 
    k_l_list = np.array([[0,0.000],[1,0.021],[2,-0.303],[3,-0.194],[4,-0.132],[5,-0.104],[6,-0.089],[7,-0.081],
                    [8,-0.076],[9,-0.072],[10,-0.069],[12,-0.064],[15,-0.058],[20,-0.051],[30,-0.040],
                    [40,-0.033],[50,-0.027],[70,-0.020],[100,-0.014],[150,-0.010],[200,-0.007]])

    f = interp1d(k_l_list[:,0],k_l_list[:,1])
    k_l = f(l)

    return k_l
