import numpy as np
import lmfit
from lmfit import minimizer, parameters, report_fit



# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, my_lines, setpoint_offset, plot_fit = false):
    a = params['a']
    a0 = params['ab0']
    a1 = params['ab1']
    a2 = params['ab2']
    w0 = params['w0']
    w1 = params['w1']
    x_offset = params['x_offset']
    y_offset = params['y_offset']
    cnt0 = params['cnt0']
    cnt1 = params['cnt1']
    cnt2 = params['cnt2']
    slope = params['slope']
   
    ampl = []
    freqs = []
    for k in range(len(my_lines)):
        ampl.append(params['a' + str(k)])

        freqs.append(my_lines[k][0]/1e6 - setpoint_offset*1e12/1e6 - x_offset)
                      
    if plot_fit == false:
        x_fit = x
    else:
        x_plot = np.linspace(np.min(x), np.max(x), 1000)
        x_fit = x_plot

    model = y_offset + slope * x_fit
    model += -a0 * np.exp( -(x_fit - cnt0)**2/(2.0*w0**2) )
    model += -a1 * np.exp( -(x_fit - cnt1)**2/(2.0*w0**2) )
    model += -a2 * np.exp( -(x_fit - cnt2)**2/(2.0*w0**2) )
    for k in range(len(my_lines)):            
        model += a * ampl[k] * np.exp( -(x_fit - freqs[k])**2/(2.0*w1**2) )
        #model += a * ampl[k] * w1**2 / ( (x_fit - freqs[k])**2 + (w1**2) )

    if plot_fit == false:
        return model - data
    else:
        return (x_plot, model)



def fit_rb(x, y, my_lines, setpoint_offset):
        params = parameters()
 
        params.add('a', value=+1.0, min=0.0, max=2.0, vary = false)
        params.add('w0', value=295.0, min=1.0, max=2000, vary = true)
        params.add('w1', value=10.48, min=1.0, max=50, vary = true)
        params.add('x_offset', value=+20.0, min=-30.0, max = 30.0, vary = true)
        params.add('y_offset', value=0.0, min=-2.0, max=2.0, vary = true)
        params.add('slope', value=-1.0/6000.0, min=-2.0, max=2.0, vary = true)
        
        params.add('ab0', value=+0.43, min=0.0, max=5.0, vary = true)
        params.add('ab1', value=+0.72, min=0.0, max=5.0, vary = true)
        params.add('ab2', value=+0.72, min=0.0, max=5.0, vary = true)
        params.add('cnt0', value=-1914.0, min=-3000.0, max=-1000, vary = true)
        params.add('cnt1', value=-785.0, min=-1000.0, max=0, vary = true)
        params.add('cnt2', value=+2000.0, min=1000.0, max=3000.0, vary = true)

         
        for k in range(len(my_lines)):
            params.add('a' + str(k), value = 0.1, min = 0.0, max = 10.0, vary = true)


        # do fit, here with leastsq model
        minner = minimizer(fcn2min, params, fcn_args=(x, y, my_lines, setpoint_offset))
        result = minner.minimize()
        
        # store the confidence data from the fit
        con_report = lmfit.fit_report(result.params)
        
        (x_plot, model) = fcn2min(result.params, x, y, my_lines, setpoint_offset, plot_fit = true)

        return (x_plot, model, result)

