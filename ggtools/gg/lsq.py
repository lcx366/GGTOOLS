import numpy as np
from scipy.optimize import curve_fit

from .fit_func import func
    
def lsqm(xdata,ydata):
    '''
    Linearly fit f(x) = a0 + a1/T*x using the Least Square algorithm.

    Usage: 
    a1,sigma_a1,a0,sigma_a0 = LSM(xdata,ydata)

    Inputs: 
    xdata -> [float array] x data
    ydata -> [float array] y data

    Outputs:
    a1 -> [float] Slope
    sigma_a1 -> [float] standard deviation for slope
    a0 -> [float] Intercept
    sigma_a0 -> [float] standard deviation for Intercept
    '''
    # estimate parameters and its covariance matrix
    if not np.diff(ydata).all():
        return 0,0,ydata[0],0

    popt,pcov = curve_fit(func, xdata, ydata)
    # calculate standard deviation
    perr = np.sqrt(np.diag(pcov))

    # values of parameters
    a0 = popt[0]
    a1 = popt[1]

    # standard deviations of parameters
    sigma_a0 = perr[0]
    sigma_a1 = perr[1]

    return a1,sigma_a1,a0,sigma_a0

def wlsqm(xdata,ydata,sigma):
    '''
    Linearly fit f(x) = a0 + a1/T*x using the Weighted Least Square algorithm.

    Usage: 
    a1,sigma_a1,a0,sigma_a0 = WLSM(xdata,ydata)

    Inputs: 
    xdata -> [float array] x data
    ydata -> [float array] y data
    sigma -> [float array] standard deviations of y data

    Outputs:
    a1 -> [float] Slope
    sigma_a1 -> [float] standard deviation for slope
    a0 -> [float] Intercept
    sigma_a0 -> [float] standard deviation for Intercept
    '''
    if not np.diff(ydata).all():
        return 0,0,ydata[0],0
    
    p0 = [1,1] # Initial guess for the parameters
    absolute_sigma = True
    
    # estimate parameters and its covariance matrix
    popt,pcov = curve_fit(func,xdata,ydata,p0,sigma,absolute_sigma)
    # calculate standard deviation
    perr = np.sqrt(np.diag(pcov))

    # values of parameters
    a0 = popt[0]
    a1 = popt[1]

    # standard deviations of parameters
    sigma_a0 = perr[0]
    sigma_a1 = perr[1]

    return a1,sigma_a1,a0,sigma_a0

def ilsqm(xdata,ydata):
    '''
    Fit f(x) = a0 + a1/T*x using the Iterative Least Square method. The 3-sigma rule is used to eliminate outliers.

    Usage: 
    a1,sigma_a1,less_than_3s,a0,sigma_a0 = ILSM(xdata,ydata)

    Inputs: 
    xdata -> [float array] x data
    ydata -> [float array] y data

    Outputs:
    a1 -> [float] Slope
    sigma_a1 -> [float] standard deviation for slope
    less_than_3s [bool array] identify whether a data point is normal(True) or abnormal(False) based on the 3-sigma rule
    a0 -> [float] Intercept
    sigma_a0 -> [float] standard deviation for Intercept
    '''
    if not np.diff(ydata).all():
        return 0,0,np.ones_like(ydata,dtype=bool),ydata[0],0

    # set a factor to multiply a std
    factor = 3.0
    n,m = len(xdata),2

    # fit the function using the xdata and ydata
    popt,pcov = curve_fit(func, xdata, ydata)
    perr = np.sqrt(np.diag(pcov))

    a0 = popt[0]
    a1 = popt[1]
    
    sigma_a0 = perr[0]
    sigma_a1 = perr[1]

    ydata_fit = func(xdata,a0,a1)
    residual = np.abs(ydata-ydata_fit) 
    rms = np.sqrt(np.dot(residual,residual)/(n-m))
    less_than_3s = residual < factor*rms

    # select data that is within 3 sigma
    ydata_i = ydata[less_than_3s]
    xdata_i = xdata[less_than_3s]

    popt_i,pcov_i = curve_fit(func, xdata_i, ydata_i) 
    perr_i = np.sqrt(np.diag(pcov_i))

    a0_i = popt_i[0] 
    a1_i = popt_i[1] 

    sigma_a0_i = perr_i[0] 
    sigma_a1_i = perr_i[1]

    ydata_fit_i = func(xdata_i,a0_i,a1_i)
    residual_i = np.abs(ydata_i-ydata_fit_i)
    n = len(residual_i)
    rms_i = np.sqrt(np.dot(residual_i,residual_i)/(n-m))
    ydata_fit_ii = func(xdata,a0_i,a1_i)
    residual_ii = np.abs(ydata-ydata_fit_ii)
    less_than_3s = residual_ii < factor*rms_i

    ydata_ii = ydata[less_than_3s] 
    xdata_ii = xdata[less_than_3s]

    while(list(xdata_i) != list(xdata_ii)):

        xdata_i = xdata_ii
        ydata_i = ydata_ii

        popt_i,pcov_i = curve_fit(func, xdata_i, ydata_i)
        perr_i = np.sqrt(np.diag(pcov_i))

        a0_i = popt_i[0]
        a1_i = popt_i[1]

        sigma_a0_i = perr_i[0]
        sigma_a1_i = perr_i[1]


        ydata_fit_i = func(xdata_i,a0_i,a1_i)
        residual_i = np.abs(ydata_i-ydata_fit_i)
        n = len(residual_i)
        rms_i = np.sqrt(np.dot(residual_i,residual_i)/(n-m))
        ydata_fit_ii = func(xdata,a0_i,a1_i)
        residual_ii = np.abs(ydata-ydata_fit_ii)
        less_than_3s = residual_ii < factor*rms_i

        ydata_ii = ydata[less_than_3s] 
        xdata_ii = xdata[less_than_3s]

    # the final parameters, variance and data
    a0 = a0_i
    a1 = a1_i

    sigma_a0 = sigma_a0_i
    sigma_a1 = sigma_a1_i
    
    factor_sigma = factor*rms_i

    return a1,sigma_a1,less_than_3s,a0,sigma_a0

def iwlsqm(xdata,ydata,sigma):
    '''
    Linearly fit f(x) = a0 + a1/T*x using the Iterative Weighted Least Square algorithm.

    Usage: 
    a1,sigma_a1,less_than_3s,a0,sigma_a0 = IWLSM(xdata,ydata,sigma)

    Inputs: 
    xdata -> [float array] x data
    ydata -> [float array] y data
    sigma -> [float array] standard deviations of y data
    
    Outputs:
    a1 -> [float] Slope
    sigma_a1 -> [float] standard deviation for slope
    less_than_3s [bool array] identify whether a data point is normal(True) or abnormal(False) based on the 3-sigma rule
    a0 -> [float] Intercept
    sigma_a0 -> [float] standard deviation for Intercept
    '''
    if not np.diff(ydata).all():
        return 0,0,np.ones_like(ydata,dtype=bool),ydata[0],0
    
    p0 = [1,1] # Initial guess for the parameters
    factor = 3 # set a factor to multiply a std
    n,m = len(xdata),len(p0)
    absolute_sigma = True
    
    # fit the function using the xdata and ydata
    popt,pcov = curve_fit(func, xdata, ydata,p0, sigma,absolute_sigma)
    perr = np.sqrt(np.diag(pcov))

    a0 = popt[0]
    a1 = popt[1]
    
    sigma_a0 = perr[0]
    sigma_a1 = perr[1]

    ydata_fit = func(xdata,a0,a1)
    residual = np.abs(ydata-ydata_fit) 
    W = 1/sigma**2
    rms = np.sqrt(np.dot(residual*W,residual)/((n-m)/n*W.sum()))
    less_than_3s = residual < factor*rms

    # select data that is in 3 sigma
    ydata_i = ydata[less_than_3s]
    xdata_i = xdata[less_than_3s]
    sigma_i = sigma[less_than_3s]
    
    popt_i,pcov_i = curve_fit(func, xdata_i, ydata_i,p0,sigma_i,absolute_sigma) 
    perr_i = np.sqrt(np.diag(pcov_i))

    a0_i = popt_i[0] 
    a1_i = popt_i[1] 

    sigma_a0_i = perr_i[0] 
    sigma_a1_i = perr_i[1]

    ydata_fit_i = func(xdata_i,a0_i,a1_i)
    residual_i = np.abs(ydata_i-ydata_fit_i)
    W_i = 1/sigma_i**2
    n = len(xdata_i)
    rms_i = np.sqrt(np.dot(residual_i*W_i,residual_i)/((n-m)/n*W_i.sum()))
    ydata_fit_ii = func(xdata,a0_i,a1_i)
    residual_ii = np.abs(ydata-ydata_fit_ii)
    less_than_3s = residual_ii < factor*rms_i

    ydata_ii = ydata[less_than_3s] 
    xdata_ii = xdata[less_than_3s]
    sigma_ii = sigma[less_than_3s]
    
    while(list(xdata_i) != list(xdata_ii)):

        xdata_i = xdata_ii
        ydata_i = ydata_ii
        sigma_i = sigma_ii
        
        popt_i,pcov_i = curve_fit(func, xdata_i, ydata_i,p0,sigma_i,absolute_sigma)
        perr_i = np.sqrt(np.diag(pcov_i))

        a0_i = popt_i[0]
        a1_i = popt_i[1]

        sigma_a0_i = perr_i[0]
        sigma_a1_i = perr_i[1]

        ydata_fit_i = func(xdata_i,a0_i,a1_i)
        residual_i = np.abs(ydata_i-ydata_fit_i)
        W_i = 1/sigma_i**2
        n = len(xdata_i)
        rms_i = np.sqrt(np.dot(residual_i*W_i,residual_i)/((n-m)/n*W_i.sum()))
        ydata_fit_ii = func(xdata,a0_i,a1_i)
        residual_ii = np.abs(ydata-ydata_fit_ii)
        less_than_3s = residual_ii < factor*rms_i

        ydata_ii = ydata[less_than_3s] 
        xdata_ii = xdata[less_than_3s]
        sigma_ii = sigma[less_than_3s]

    # the final parameters, variance and data
    a0 = a0_i
    a1 = a1_i

    sigma_a0 = sigma_a0_i
    sigma_a1 = sigma_a1_i
    
    factor_sigma = factor*rms_i

    return a1,sigma_a1,less_than_3s,a0,sigma_a0