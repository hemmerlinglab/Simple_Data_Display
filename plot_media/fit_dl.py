import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit

import matplotlib.pyplot as plt
from fit_dunham import energy

c = 299792458
kB = 1.3806482e-23
h_planck = 6.63607e-34
amu = 1.66e-27

def get_doppler(T, f0):

    return np.sqrt(8*kB*T*np.log(2)/((35+27)*amu*c**2)) * f0


def gauss(x, x0, A, w):

    return A * np.exp( -(x-x0)**2/(2*w**2) )


# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, plot_fit = False):
    a1 = params['a1']
    a2 = params['a2']
    x1 = params['x1']
    x2 = params['x2']
    w1 = params['w1']
    w2 = params['w2']
    T1 = params['T1']
    T2 = params['T2']
   

    cnt_freq = 3*382.11035e12/100.0/c

    w1 = get_doppler(T1, 100*c*cnt_freq)
    w2 = get_doppler(T2, 100*c*cnt_freq)

    nus = x*1e6 # in Hz

    df = 1e9

    if plot_fit == False:
        xfit = nus
    else:
        xfit = np.linspace(np.min(nus) - df, np.max(nus) + df, 500)

    y_th = gauss(xfit, x1, a1, w1) + a2 * gauss(xfit, x2, a2, w2)


    if plot_fit == False:
        return (y_th - data)
    else:
        return (xfit/1e6, y_th)


def fit_dl(x, y):
        params = Parameters()
 
        params.add('a1', value=1.5, min=0.0, max=5.0, vary = True)
        params.add('a2', value=1.5, min=0.0, max=5.0, vary = True)
        params.add('w1', value=0.1, min=0.0, max=1.0e9, vary = True)
        params.add('w2', value=0.1, min=0.0, max=1.0e9, vary = True)
        params.add('x1', value=0.1, min=-10e9, max=10.0e9, vary = True)
        params.add('x2', value=1.5e9, min=-10e9, max=10.0e9, vary = True)
        params.add('T1', value=20.3, min=0.1, max=50.0, vary = True)
        params.add('T2', value=20.3, min=0.1, max=50.0, vary = True)
         
        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(x, y))
        result = minner.minimize()
        
        # Store the Confidence data from the fit
        #con_report = lmfit.fit_report(result.params)
        
        (x_fit, y_fit) = fcn2min(result.params, x, y, plot_fit = True)

        #return (x_plot, model, result)

        return (x_fit, y_fit, result)

