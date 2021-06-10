"""
Utility functions.

"""

import numpy as np

__author__ = "Nong Artrith, Alexander Urban"
__email__ = "nartrith@atomistic.net"
__maintainer__ = "Nong Artrith, Alexander Urban"
__maintainer_email__ = ("nartrith@atomistic.net"
                        + ", aurban@atomistic.net")
__date__ = "2021-06-08"
__version__ = "0.1"


def gauss_function(x, x0, s):
    """
    Normalized Gaussian

    Args:
      x (float): argument
      x0 (float): mean
      s (float): standard deviation

    """
    return np.exp(-0.5*((x-x0)/s)**2)/(s*np.sqrt(2.0*np.pi))


def lorentz_function(x, x0, w):
    """
    Normalized Lorentz (Cauchy) function

    Args:
      x (float): argument
      x0 (float): center
      w (float): full width at half maximum (FWHM)

    """
    return 2.0/(w*np.pi)/(1.0 + 4.0*(x - x0)**2/w**2)


def broaden(x_in, y_in, fwhm, fwhm2=None, xlim=None, n=100, sigma=False,
            lorentz=False):
    """
    Broaden a line shape via convolution with a distribution function.
    A Gaussian function is used unless `lorentz` is set to True, in
    which case a Lorentzian function is used.

    Args:
      x_in, y_in (ndarray): X and Y values of the line shape
      fwhm (float): full width at half maximum (FWHM) of the distribution
        function
      fwhm2 (float): if specified, the FWHM is linearly varied from
        `fwhm` to `fwhm2` from the lowest to the highest X value
      xlim (tuple): Tuple (x0, x1) with the minimal and maximal X values
      n (int): number of discretization points
      sigma (bool): if True, assume that `fwhm` is the standard deviation
      lorentz (bool): if True, perform convolution with a Lorentzian function
        instead of a Gaussian function

    Returns:
      x, y (ndarray)

    """

    if fwhm2 is None:
        scale = np.ones(n)
    else:
        scale = np.linspace(1.0, fwhm2/fwhm, n)

    if lorentz:
        func = lorentz_function
        sigma = True
    else:
        func = gauss_function

    if sigma or lorentz:
        s = fwhm
    else:
        s = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))

    x0 = np.min(x_in)
    x1 = np.max(x_in)
    if xlim is not None:
        if xlim[0] is not None:
            x0 = max(xlim[0], x0)
        if xlim[1] is not None:
            x1 = min(xlim[1], x1)

    # determine non-zero y values so that we can ignore exact zero
    # values at beginning and end of x range
    idx = [i for i in range(len(x_in)) if y_in[i] != 0.0]

    y = np.zeros(n)
    x = np.linspace(x0, x1, n)
    y[0] = y_in[0]
    y[-1] = y_in[-1]
    for i in range(max(idx[0], 1), min(idx[-1], len(x_in)-1)):
        dx = 0.5*(x_in[i+1]-x_in[i-1])
        for j in range(1, n-1):
            y[j] += y_in[i]*func(x[j], x_in[i], scale[j]*s)*dx

    return x, y
