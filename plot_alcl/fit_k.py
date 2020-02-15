import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit

# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, plot_fit = False):
    a = params['a']
    w = params['w']
    x_offset = params['x_offset']
    y_offset = params['y_offset']
    
    iso_abund = 1.0/100.0 * np.array([93.2581, 93.2581, 6.7302, 6.7302])
    
    freqs = 173.1-14.4+np.array([+288.6-16.1, -173.1+14.4, +236.2+158.8-8.4, +236.2-95.3+8.4])
    
    ampl = []
    for k in range(len(iso_abund)):
        ampl.append(params['a' + str(k)])

    if plot_fit == False:
        freqs = freqs - x_offset

        model = y_offset

        for k in range(len(iso_abund)):
            model += a * ampl[k] * iso_abund[k] * np.exp( -(x - freqs[k])**2/(2.0*w**2) )
 
        return model - data
    else:
        x_plot = np.linspace(np.min(x), np.max(x), 200)
        
        freqs = freqs - x_offset

        model = y_offset
        for k in range(len(iso_abund)):
            model += a * ampl[k] * iso_abund[k] * np.exp( -(x_plot - freqs[k])**2/(2.0*w**2) )
        
        return (x_plot, model)


def fit_k(x, y):
        params = Parameters()
 
        params.add('a', value=-5.0, min=-10.0, max=0.0, vary = True)
        params.add('w', value=50.0, min=1.0, max=2000, vary = True)
        params.add('x_offset', value=50, min=np.min(x), max = np.max(x), vary = True)
        params.add('y_offset', value=0.0, min=-2.0, max=2.0, vary = True)

         
        iso_abund = np.array([93.2581, 93.2581, 6.7302, 6.7302])
        for k in range(len(iso_abund)):
            params.add('a' + str(k), value = 1.0, min = 0.0, max = 10.0, vary = True)


        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(x, y))
        result = minner.minimize()
        
        # Store the Confidence data from the fit
        con_report = lmfit.fit_report(result.params)
        
        (x_plot, model) = fcn2min(result.params, x, y, plot_fit = True)

        return (x_plot, model, result)

