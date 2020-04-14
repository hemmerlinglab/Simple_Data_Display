import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit

from helper import get_spec_lines, scale_coeff, get_reduced_dunham, get_dunham

def abs_sig(x, p):

    return p['p0'] + p['p1'] * np.exp( -(x - p['p2']) / p['p3'] ) * (x - p['p2'])

def gauss(x, p):

    return p['p0'] + p['p1'] * np.exp( -(x - p['p2'])**2/(2.0*p['p3']**2) )

def multi_gauss(x, p):

    f = p['p0']

    for n in range(4):
        f += p['p1'+str(n)] * np.exp( -(x - p['p2'+str(n)])**2/(2.0*p['p3']**2) )

    return f



def P_lines(x, p):

    Jmax = 10
    T = 5 # Kelvin

    (Ug, Ue) = get_reduced_dunham()
    
    # we vary the mass reduced coefficients
    Ue[0][0] = p['p0'] # TeA
    Ue[0][1] = p['p1'] # BeA
    Ue[1][1] = p['p2'] # 
    Ue[0][2] = p['p3'] # De 
    Ue[1][0] = p['p4'] # we 

    (Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)

    # vibrational states
    ve = 0
    vg = 0
    (nus35, spec35) = get_spec_lines(Yg35, Ye35, ve, vg, x, line_type = 'P', T = T, Jmax = Jmax, real_amplitudes = False)
    (nus37, spec37) = get_spec_lines(Yg37, Ye37, ve, vg, x, line_type = 'P', T = T, Jmax = Jmax, real_amplitudes = False)

    f = spec35 + spec37

    return f


 


def print_results(m):

    for k in m.params.keys():

        print('Parameter ' + k + ' = ' + str(m.params[k].value))

def fcn2min(params, x, data, func = None, plot_fit = False):

    if plot_fit == False:
        x_fit = x
    else:
        x_fit = np.linspace(np.min(x), np.max(x), 500)

    model = func(x_fit, params)
    #y_offset + a * np.exp( -(x_fit - x0)**2/(2.0*w**2) )
    
    if plot_fit == False:
        return model - data
    else:
        return (x_fit, model)


def fit_func(x, y, params, func):
                 
        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(x, y, func))
        result = minner.minimize()
        
        ## Store the Confidence data from the fit
        #con_report = lmfit.fit_report(result.params)
        
        (x_plot, model) = fcn2min(result.params, x, y, func = func, plot_fit = True)

        return (x_plot, model, result)


if __name__ == '__main__':

    params = Parameters()
 
    params.add('p0', value=2.4, min=-2.0, max=4.0, vary = True)
    params.add('p1', value=100.0, min=0.0, max=2000.0, vary = True)
    params.add('p2', value=0.0, min=-2, max=2, vary = True)
    params.add('p3', value=2.3, min=-2.0, max=4.0, vary = True)
         
    x = np.linspace(-10,10,100)

    y = 2.4 + 0.2 * np.exp(-x**2/(2*2.3**2))

    (x_fit, y_fit, result) = fit_func(x, y, params, gauss)

    print(result.params)

    import matplotlib.pyplot as plt
    plt.plot(x, y, 'o')
    plt.plot(x_fit, y_fit)
    plt.show()
