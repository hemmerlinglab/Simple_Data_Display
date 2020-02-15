import lmfit
import numpy as np
from lmfit import Minimizer, Parameters, report_fit

amu = 1.660539e-27
c = 299792458
kB = 1.38064852e-23

# define objective function: returns the array to be minimized
def fcn2min(params, x, data, return_plot = False):
    """Model a decaying sine wave and subtract data."""
    a1 = params['a1']
    a2 = params['a2']
    a3 = params['a3']
    a4 = params['a4']
    a5 = params['a5']
    a6 = params['a6']
    a7 = params['a7']
    a8 = params['a8']
    a9 = params['a9']
    a10 = params['a10']
    a11 = params['a11']
    a12 = params['a12']
    shift = params['shift']
    T = params['T']
    
    m39 = 39 * amu
    m41 = 41 * amu
    
    if return_plot:
        x = np.linspace(np.min(x), np.max(x), 250)

    k39d2line = 391.01617e12

    k39_23 = k39d2line - 173.1e6 + 14.4e6
    k39_12 = k39d2line + 288.6e6 - 6.7e6
    k39_11 = k39d2line + 288.6e6 - 16.1e6
    k39_10 = k39d2line + 288.6e6 - 19.4e6
    k39_22 = k39d2line - 173.1e6 - 6.7e6
    k39_21 = k39d2line - 173.1e6 - 16.1e6
    
    k41_12 = k39d2line + 236.2e6 + 158.8e6 - 5.0e6
    k41_11 = k39d2line + 236.2e6 + 158.8e6 - 8.4e6
    k41_10 = k39d2line + 236.2e6 + 158.8e6 - 8.4e6
    k41_21 = k39d2line + 236.2e6 - 95.3e6 - 8.4e6
    k41_22 = k39d2line + 236.2e6 - 95.3e6 - 5.0e6
    k41_23 = k39d2line + 236.2e6 - 95.3e6 + 8.4e6

    model  = a1 * (c/(k39_23*np.sqrt(2*kB*T/m39))) * np.exp( -m39 * c**2 * (x + shift - k39_23)**2/(2*kB*T* k39_23**2) )
    model += a2 * (c/(k39_12*np.sqrt(2*kB*T/m39))) * np.exp( -m39 * c**2 * (x + shift - k39_12)**2/(2*kB*T* k39_12**2) )
    model += a3 * (c/(k39_11*np.sqrt(2*kB*T/m39))) * np.exp( -m39 * c**2 * (x + shift - k39_11)**2/(2*kB*T* k39_11**2) )
    model += a4 * (c/(k39_10*np.sqrt(2*kB*T/m39))) * np.exp( -m39 * c**2 * (x + shift - k39_10)**2/(2*kB*T* k39_10**2) )
    model += a5 * (c/(k41_12*np.sqrt(2*kB*T/m41))) * np.exp( -m41 * c**2 * (x + shift - k41_12)**2/(2*kB*T* k41_12**2) )
    model += a6 * (c/(k39_21*np.sqrt(2*kB*T/m39))) * np.exp( -m39 * c**2 * (x + shift - k39_21)**2/(2*kB*T* k39_21**2) )
    model += a7 * (c/(k39_22*np.sqrt(2*kB*T/m39))) * np.exp( -m39 * c**2 * (x + shift - k39_22)**2/(2*kB*T* k39_22**2) )
    model += a8 * (c/(k41_11*np.sqrt(2*kB*T/m41))) * np.exp( -m41 * c**2 * (x + shift - k41_11)**2/(2*kB*T* k41_11**2) )
    model += a9 * (c/(k41_10*np.sqrt(2*kB*T/m41))) * np.exp( -m41 * c**2 * (x + shift - k41_10)**2/(2*kB*T* k41_10**2) )
    model += a10 * (c/(k41_21*np.sqrt(2*kB*T/m41))) * np.exp( -m41 * c**2 * (x + shift - k41_21)**2/(2*kB*T* k41_21**2) )
    model += a11 * (c/(k41_22*np.sqrt(2*kB*T/m41))) * np.exp( -m41 * c**2 * (x + shift - k41_22)**2/(2*kB*T* k41_22**2) )
    model += a12 * (c/(k41_23*np.sqrt(2*kB*T/m41))) * np.exp( -m41 * c**2 * (x + shift - k41_23)**2/(2*kB*T* k41_23**2) )

    if return_plot == False:
        return model - data
    else:
        return x, model


def my_fit(x, y):


    # create a set of Parameters
    params = Parameters()
    params.add('a1', value=1.0, min = 0)
    params.add('a2', value=1.0, min = 0)
    params.add('a3', value=1.0, min = 0)
    params.add('a4', value=1.0, min = 0)
    params.add('a5', value=1.0, min = 0)
    params.add('a6', value=1.0, min = 0)
    params.add('a7', value=1.0, min = 0)
    params.add('a8', value=1.0, min = 0)
    params.add('a9', value=1.0, min = 0)
    params.add('a10', value=1.0, min = 0)
    params.add('a11', value=1.0, min = 0)
    params.add('a12', value=1.0, min = 0)
    params.add('shift', value=-5.0e6, min = -500.0e6, max = 500.0e6)
    params.add('T', value=10.0, min = 0.1, max = 70.0)

    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(x, y))
    result = minner.minimize()

    return result



