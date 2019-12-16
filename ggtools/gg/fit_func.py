def func(x,a0,a1):
    '''
    Define a linear function f(x) = a0 + a1/T*x to be fitted, where a0 and a1 are parameters for intercept and slope.

    Usage:
    y = func(x,a0,a1)

    Inputs:
    x -> [float array] Independent variables
    a0 -> [float] Intercept
    a1 -> [float] Slope

    Outputs:
    y -> [float array] Dependent variables

    Examples:
    >>> import numpy as np
    >>> xs = np.arange(0,120,12)
    >>> a0,a1 = 1,2
    >>> ys = func(xs,a0,a1)
    >>> print(xs)
    [  0  12  24  36  48  60  72  84  96 108]
    >>> print(ys)
    [ 1.  3.  5.  7.  9. 11. 13. 15. 17. 19.]
    '''
    T = 12 # 12 months in one year
    return a0 + a1/T*x