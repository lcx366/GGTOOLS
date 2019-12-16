import numpy as np
from numpy.linalg import inv,norm
import matplotlib.pyplot as plt

def solve_lambda_x(A,y,lambs_exponent=np.arange(-10,5,0.01)):
    '''
    Minimize the Lagrangian L(x,lambda) = || A*x - y||**2 + lambda*||x||**2 based on the idea of Tikhonov Regularization.
    This program is used to roughly estimate the parameters x and the Lagrange multiplier lambda. 
    The L-curve method is applied to get the proper Lagrangian multiplier.

    Usage:
    estimate_lamb,estimate_x,log10_residual_norm,log10_solution_norm,curvature = solve_lambda_x(A,y) 

    Inputs:
    A -> [float 2d array] Design matrix
    y -> [float array] Measurements

    Parameters:
    lambs_exponent -> [optional, float 3d/4d array, default = np.arange(-10,5,0.01)] Exponent for lambda with base of 10
    
    Outputs:
    estimate_lamb -> [float] Lagrange multiplier
    estimate_x -> [float array] Estimated parameters
    log10_residual_norm -> [float array] log10(||A*x-y||) with lambda taking 10**lambs_exponent
    log10_solution_norm -> [float array] log10(||x||) with lambda taking 10**lambs_exponent
    curvature -> [float array] curvature of the L-curve, where the ordinate of the curve is log10_solution_norm and the abscissa is log10_residual_norm.

    For more information, please refer to 
    (1) [NumPy/SciPy Recipes for Data Science: Regularized Least Squares Optimization](https://www.researchgate.net/publication/274138835_NumPy_SciPy_Recipes_for_Data_Science_Regularized_Least_Squares_Optimization)
    (2) [Choosing the Regularization Parameter](http://www2.compute.dtu.dk/~pcha/DIP/chap5.pdf)
    '''  
    np.seterr(divide='ignore',invalid='ignore')
    m = A.shape[1]
  
    # set a series of Lagrangian multiplier       
    log10_residual_norm,log10_solution_norm = [],[]
    lambs = np.float_power(10,lambs_exponent)
    for lamb in lambs:
        x = np.dot(inv(np.dot(A.T,A)+lamb*np.eye(m)),np.dot(A.T,y))
        residual_norm = norm(np.dot(A,x)-y)
        solution_norm = norm(x)
        log10_residual_norm.append(np.log10(residual_norm))
        log10_solution_norm.append(np.log10(solution_norm))
    log10_residual_norm = np.array(log10_residual_norm)
    log10_solution_norm = np.array(log10_solution_norm)  
    
    # calculate the curvature of the L-curve
    g1 = np.gradient(log10_solution_norm, log10_residual_norm)
    g1[np.isnan(g1)] = -np.inf
    g2 = np.gradient(g1,log10_residual_norm)
    g2[np.isnan(g2)] = np.inf
    curvature = np.abs(g2)/(1+g1**2)**1.5
    curvature[np.isnan(curvature)] = 0
    curvature[np.isinf(curvature)] = 0
    index_curvature_max = np.argmax(curvature)
    estimate_lamb = lambs[index_curvature_max]
    estimate_x = np.dot(inv(np.dot(A.T,A)+estimate_lamb*np.eye(m)),np.dot(A.T,y))
    return estimate_lamb,estimate_x,log10_residual_norm,log10_solution_norm,curvature

def L_curve(A,y,visible=None):
    '''
    Minimize the Lagrangian L(x,lambda) = || A*x - y||**2 + lambda*||x||**2 based on the idea of Tikhonov Regularization.
    This program is used to accurately estimate the parameters x and the Lagrange multiplier lambda. 
    The final Lagrange multiplier is determained by the L-curve method. The L-curve can be visualized by outputing an image.

    Usage:
    accu_lamb,accu_x = L_curve(A,y) 

    Inputs:
    A -> [float 2d array] Design matrix
    y -> [float array] Measurements

    Parameters:
    visible -> [optional, str, default = None] If None, the visualization of L-vurve will be closed. If 'visible', the L-curve will be visualized by outputing an image.
    
    Outputs:
    accu_lamb -> [float] Lagrange multiplier
    accu_x -> [float array] Estimated parameters
    
    For more information, please refer to 
    (1) [NumPy/SciPy Recipes for Data Science: Regularized Least Squares Optimization](https://www.researchgate.net/publication/274138835_NumPy_SciPy_Recipes_for_Data_Science_Regularized_Least_Squares_Optimization)
    (2) [Choosing the Regularization Parameter](http://www2.compute.dtu.dk/~pcha/DIP/chap5.pdf)
    '''  
    m = A.shape[1]
            
    # Estimate the Lagrange multiplier roughly  
    appr_lamb,appr_x,log10_residual_norm,log10_solution_norm,appr_curvature = solve_lambda_x(A,y)
    
    # Estimate the Lagrange multiplier accurately
    lambs_exponent = np.linspace(np.log10(appr_lamb)-2,np.log10(appr_lamb)+2,2000)
    accu_lamb,accu_x,log10_residual_norm,log10_solution_norm,accu_curvature = solve_lambda_x(A,y,lambs_exponent)

    if visible is not None:
        fig_dir = 'figures/'
        if not os.path.exists(fig_dir): os.makedirs(fig_dir) 
        # plot
        plt.clf()
        fig, (ax1, ax2) = plt.subplots(1, 2,dpi=200)
        # make a little extra space between the subplots
        fig.subplots_adjust(wspace=0.4)
        ax1.plot(log10_residual_norm,log10_solution_norm)
        ax1.set_xlabel(r'$\log \parallel A x_{\lambda}-y \parallel_2$')
        ax1.set_ylabel(r'$\log \parallel x_{\lambda} \parallel_2$')
        ax1.set_title('L-Curve')
        ax2.plot(lambs_exponent,accu_curvature)
        ax2.set_xlabel(r'$\log \parallel \lambda \parallel_2$')
        ax2.set_ylabel('Curvature')
        ax2.set_title('curvature of L-Curve')
        plt.savefig(fig_dir+'L-Curve.png')
    return accu_lamb,accu_x