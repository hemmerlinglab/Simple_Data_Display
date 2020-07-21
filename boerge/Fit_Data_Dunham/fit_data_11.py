import numpy as np
import matplotlib.pyplot as plt
import pickle
import lmfit
from lmfit import Minimizer, Parameters, report_fit

from aux_functions import *

from get_spectrum import get_spectrum

###########################################################################
# this program fits the transitions and extracts the Dunham coefficients
###########################################################################

def get_data(filename):
    arr = pickle.load( open( filename, "rb" ) )

    x = arr['f_arr']
    y = arr['sig_arr']
    line_info = arr['data']

    return (x, y, line_info)

def fcn2min(params, x, data, func = None, plot_fit = False):

    y_offset = params['p0']
    a = params['p1']
    x0 = params['p2']
    w = params['p3']

    if plot_fit == False:
        x_fit = x
    else:
        x_fit = np.linspace(np.min(x), np.max(x), 500)

    model = y_offset + a * np.exp( -(x_fit - x0)**2/(2.0*w**2) )
    
    if plot_fit == False:
        return model - data
    else:
        return (x_fit, model)


def fit_line(x, y):

    y_min = np.min(y)
    y_max = np.max(y)
    x_mean = np.mean(x)
    x_min = np.min(x)
    x_max = np.max(x)

    params = Parameters()

    params.add('p0', value=y_min, min=-1.0, max=1.0, vary = True)
    params.add('p1', value=(y_max - y_min), min=0.0, max=3.0, vary = True)
    params.add('p2', value=x_mean, min=x_mean-10e9, max=x_mean+10e9, vary = True)
    params.add('p3', value=(x_max-x_min)/10.0, min=100e6, max=(x_max-x_min), vary = True)

    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(x, y, None))
    result = minner.minimize()
    
    (xf, yf) = fcn2min(result.params, x, y, plot_fit = True)
    
    return (result, xf, yf)






######################################################
# main
######################################################

# Q00 line
(x_arr, y_arr, line_info) = get_data('../20200619/data.pickle')

height_arr = []
data = []
line_type = []

vg = 0
ve = 0



###########################################
# apply wavemeter offset to data

wavemeter_offset = 3 * -15.0e6

for k in range(len(x_arr)):
    x_arr[k] += wavemeter_offset

###########################################




for k in range(len(x_arr)):

    Jg = line_info[k]['Jg']
    Je = line_info[k]['Je']

    line_type.append(get_line_type(Jg, Je))
    
    x = x_arr[k]
    y = y_arr[k]

    #if not 'skip' in line_info[k].keys():
    if True:
        (result, xf, yf) = fit_line(x, y)

        plt.plot(xf, yf, '-')
        plt.plot(x, y, '.-')

        plt.text(np.mean(xf), np.max(yf), str(k))

        height_arr.append(result.params['p1'].value)

        # create data matrix to fit Dunham coefficients
        hlp = [vg, Jg, ve, Je, result.params['p2'].value, line_info[k]['iso']]

        data.append(hlp)


data2losalamos(data, save_filename = 'data_11.csv')

plt.show()




