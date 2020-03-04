import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit

# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, y1, y2):
    f = params['f']

    to_minimize = (np.mean(y1 - f*y2))**2
   
    return to_minimize



def fit_factor(x, y1, y2):
        params = Parameters()
 
        params.add('f', value=+1.0, min=-10.0, max=10.0, vary = True)

        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(x, y1, y2))
        result = minner.minimize()
        
        # store the confidence data from the fit
        con_report = lmfit.fit_report(result.params)
        
        return (result)

