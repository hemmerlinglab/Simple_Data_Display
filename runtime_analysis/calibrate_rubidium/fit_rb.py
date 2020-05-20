import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit



# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, my_lines, plot_fit = False):
    a = params['a']
    a0 = params['ab0']
    a1 = params['ab1']
    w0 = params['w0']
    w1 = params['w1']
    x_offset = params['x_offset']
    y_offset = params['y_offset']
    cnt0 = params['cnt0']
    cnt1 = params['cnt1']
   
    ampl = []
    freqs = []
    for k in range(len(my_lines)):
        ampl.append(params['a' + str(k)])

        freqs.append(my_lines[k][0])
                      
    if plot_fit == False:
        x_fit = x
    else:
        x_plot = np.linspace(np.min(x), np.max(x), 1000)
        x_fit = x_plot

    model = y_offset
    model += -a0 * np.exp( -(x_fit - cnt0)**2/(2.0*w0**2) )
    model += -a1 * np.exp( -(x_fit - cnt1)**2/(2.0*w0**2) )
    for k in range(len(my_lines)):            
        model += a * ampl[k] * np.exp( -(x_fit - freqs[k] - x_offset)**2/(2.0*w1**2) )
        #model += a * ampl[k] * w1**2 / ( (x_fit - freqs[k])**2 + (w1**2) )

    if plot_fit == False:
        return model - data
    else:
        return (x_plot, model)



def fit_rb(x, y, my_lines):
        params = Parameters()
 
        cnt_rb = 384.230115e12

        params.add('a', value=+0.3, min=0.0, max=2.0, vary = False)
        params.add('w0', value=295.0e6, min=1.0e6, max=300e6, vary = True)
        params.add('w1', value=10.48e6, min=1.0e6, max=100e6, vary = False)
        params.add('x_offset', value=-15.0e6, min=-100.0e6, max = 100.0e6, vary = False) # wavemeter offset
        params.add('y_offset', value=0.0, min=-2.0, max=2.0, vary = True)
        
        params.add('ab0', value=+0.23, min=0.0, max=2.0, vary = True)
        params.add('ab1', value=+0.92, min=0.0, max=2.0, vary = True)
        params.add('cnt0', value=cnt_rb-1914.0e6, min=cnt_rb-3000.0e6, max=cnt_rb-1000e6, vary = True)
        params.add('cnt1', value=cnt_rb-785.0e6, min=cnt_rb-1000.0e6, max=cnt_rb+0e6, vary = True)

         
        for k in range(len(my_lines)):
            params.add('a' + str(k), value = 0.2, min = 0.0, max = 5.0, vary = False)

        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(x, y, my_lines))
        result = minner.minimize()
        
        # Store the Confidence data from the fit
        con_report = lmfit.fit_report(result.params)
        
        (x_plot, model) = fcn2min(result.params, x, y, my_lines, plot_fit = True)

        return (x_plot, model, result)

