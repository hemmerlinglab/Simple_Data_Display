import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit

# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, get_fit = False):
    
    Te_1 = params['Te_1']
    Te_2 = params['Te_2']
    
    we_1 = params['we_1']
    we_2 = params['we_2']
    
    wexe_1 = params['wexe_1']
    wexe_2 = params['wexe_2']
    
    weye_1 = params['weye_1']
    weye_2 = params['weye_2']

    weze_1 = params['weze_1']
    weze_2 = params['weze_2']



    # x = [[1,2], [3,2], [1,0]]

    
    x2 = np.array(list(map(lambda t : t[0], x)))
    x1 = np.array(list(map(lambda t : t[1], x)))

    model  = Te_2 + (we_2 - wexe_2 * (x2 + 0.5) + weye_2 * (x2 + 0.5)**2 - weze_2 * (x2 + 0.5)**3 ) * (x2 + 0.5)
    model -= Te_1 + (we_1 - wexe_1 * (x1 + 0.5) + weye_1 * (x1 + 0.5)**2 - weze_1 * (x1 + 0.5)**3 ) * (x1 + 0.5)

    if get_fit == False:

        return model - data

    else:

        return model


def do_fit(x, y):
        params = Parameters()
 
        params.add('Te_1', value=0.0, min=0.0, max=340.0, vary = False)
        params.add('Te_2', value=38237.0, min=0.0, max=40000.0, vary = True)

        params.add('we_1', value=480.0, min=0.0, max=500.0, vary = True)
        params.add('we_2', value=440.0, min=0.0, max=500.0, vary = True)
 
        params.add('wexe_1', value=2.037, min=0.0, max=5.0, vary = True)
        params.add('wexe_2', value=2.81, min=0.0, max=8.0, vary = True)
    
        params.add('weye_1', value=0.5, min=0.0, max=5.0, vary = True)
        params.add('weye_2', value=0.1, min=0.0, max=8.0, vary = True)
 
        params.add('weze_1', value=0.5, min=0.0, max=5.0, vary = True)
        params.add('weze_2', value=0.1, min=0.0, max=8.0, vary = True)
      
        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(x, y))
        result = minner.minimize()
        
        # Store the Confidence data from the fit
        con_report = lmfit.fit_report(result.params)
        
        return (result)

