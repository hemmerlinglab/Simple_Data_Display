import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit

# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, plot_fit = False):
    a = params['a']
    w = params['w']
    x0 = params['x0']
    y_offset = params['y_offset']

    if plot_fit == False:
        model = y_offset + a * np.exp( -(x - x0)**2/(2.0*w**2) )
 
        return model - data
    else:
        x_plot = np.linspace(np.min(x), np.max(x), 200)
        
        model = y_offset + a * np.exp( -(x_plot - x0)**2/(2.0*w**2) )

        return (x_plot, model)


def fit_line(x, y):
        params = Parameters()
 
        params.add('a', value=-0.1, min=-10.0, max=0.0, vary = True)
        params.add('w', value=100.0, min=1.0, max=2000.0, vary = True)
        params.add('x0', value=-400.0, min=-2e3, max=2e3, vary = True)
        params.add('y_offset', value=1.0, min=-2.0, max=2.0, vary = True)
         
        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(x, y))
        result = minner.minimize()
        
        # Store the Confidence data from the fit
        con_report = lmfit.fit_report(result.params)
        
        (x_plot, model) = fcn2min(result.params, x, y, plot_fit = True)

        return (x_plot, model, result)

