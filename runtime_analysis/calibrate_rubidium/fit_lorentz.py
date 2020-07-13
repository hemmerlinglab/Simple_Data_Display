import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit


# define objective function: returns the array to be minimized
# function to minimize
def fcn2min_lorentz(params, x, data, plot_fit = False):
    a0 = params['a0']
    w0 = params['w0']
    cnt0 = params['cnt0']
    offset = params['offset']
   
    if plot_fit == False:
        x_fit = x
    else:
        no_points = 1000
        x_plot = np.linspace(np.min(x), np.max(x), no_points)
        x_fit = x_plot

    #model = offset + a0 * w0**2 / ( (x_fit - cnt0)**2 + (w0**2) )
    
    model = offset + a0 * np.exp( -(x_fit - cnt0)**2/(2*w0**2) )

    if plot_fit == False:
        return model - data
    else:
        return (x_plot, model)



def fit_lorentz(x_fit, y_fit, cnt, dn = 5):


        x = x_fit[cnt - dn: cnt + dn]
        y = y_fit[cnt - dn: cnt + dn]

        params = Parameters()
 
        x_min = np.min(x)
        x_max = np.max(x)
        y_min = np.min(y)
        y_max = np.max(y)
        x_mean = np.mean(x)
        y_mean = np.mean(y)

        params.add('w0',     value=(x_max - x_min)/10.0, min=(x_max-x_min)/20.0, max=2*(x_max-x_min), vary = True)
        params.add('offset', value=0.1, min=0.0, max = 1.0, vary = True)
        params.add('a0',     value=y_max - y_min, min=0.0, max=2.0, vary = True)
        params.add('cnt0',   value=x_fit[cnt], min=x_min, max=x_max, vary = True) # initial value equals the peak that was found

        # do fit, here with leastsq model
        minner = Minimizer(fcn2min_lorentz, params, fcn_args=(x, y))
        result = minner.minimize()
        
        # Store the Confidence data from the fit
        con_report = lmfit.fit_report(result.params)
        
        (x_plot, model) = fcn2min_lorentz(result.params, x, y, plot_fit = True)
        
        # get residuals
        (residuals) = fcn2min_lorentz(result.params, x, y)

        #print(result.params)

        #print(con_report)

        return (x_plot, model, result, residuals)

