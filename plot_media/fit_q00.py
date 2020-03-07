import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit

import matplotlib.pyplot as plt
from fit_dunham import energy

c = 299792458
kB = 1.3806482e-23
h_planck = 6.63607e-34
amu = 1.66e-27


def get_pop(J, B, T, N0 = 1):

    #pop = N0 * (2*J+1) * np.exp(-B * h_planck * c * J*(J+1)/(kB*T))
    pop = N0 * (2*J+1) * np.exp(-B * h_planck * J*(J+1)/(kB*T))

    #return pop / (np.sqrt(kB*T/(2*h_planck*c*B)) - 1/2)
    return pop / (np.sqrt(kB*T/(2*h_planck*B)) - 1.0/2.0)

def get_doppler(T, f0):

    return np.sqrt(8*kB*T*np.log(2)/((35+27)*amu*c**2)) * f0


def gauss(x, x0, A, w):

    return A * np.exp( -(x-x0)**2/(2*w**2) )


# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, Yg35, Ye35, plot_fit = False):
    a = params['a']
    Trot = params['Trot']
    Tcoll = params['Tcoll']
    x0 = params['x0']

    Ye = Ye35
    Yg = Yg35

    vg = 0
    ve = 0

    Jmax = 10
    cnt_freq = 3*382.11035e12/100.0/c

    Jg_arr = np.arange(0, Jmax+1)    

    w = get_doppler(Tcoll, 100*c*cnt_freq)

    nus = x*1e6 - x0 + 100 * c * cnt_freq # in Hz

    df = 1e9

    if plot_fit == False:
        xfit = nus
    else:
        xfit = np.linspace(np.min(nus) - df, np.max(nus) + df, 500)

    spectrum_Q = np.zeros(len(xfit))
    for Jg in Jg_arr:
        Je = Jg            
    
        eng = 100*c*(energy(Ye, ve, Je) - energy(Yg, vg, Jg))
    
        # apply population of ground states
        A = get_pop(Jg, 100 * c * np.abs(Yg[1][0]), Trot)
    
        spectrum_Q += gauss(xfit, eng, A, w)

    y_th = a * spectrum_Q # (0.76 * spectrum35  + 0.24 * spectrum37)


    if plot_fit == False:
        return y_th - data
    else:
        return ((xfit - 100*c*cnt_freq + x0)/1e6, y_th)


def fit_q00(x, y, Yg35, Ye35):
        params = Parameters()
 
        params.add('a', value=0.1, min=0.0, max=1.0, vary = True)
        params.add('Trot', value=10.3, min=0.1, max=10.0, vary = True)
        params.add('Tcoll', value=10.3, min=1e-2, max=100.0, vary = True)
        params.add('x0', value=15e9, min=-30e9, max=30e9, vary = True)
         
        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(x, y, Yg35, Ye35))
        result = minner.minimize()
        
        # Store the Confidence data from the fit
        #con_report = lmfit.fit_report(result.params)
        
        (x_fit, y_fit) = fcn2min(result.params, x, y, Yg35, Ye35, plot_fit = True)

        #return (x_plot, model, result)

        return (x_fit, y_fit, result)

