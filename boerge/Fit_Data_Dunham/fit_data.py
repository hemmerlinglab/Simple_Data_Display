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
(x_arr, y_arr, line_info) = get_data('../20200616/data.pickle')

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

    if not 'skip' in line_info[k].keys():
        (result, xf, yf) = fit_line(x, y)

        plt.plot(xf, yf, '-')
        plt.plot(x, y, '.-')

        plt.text(np.mean(xf), np.max(yf), str(k))

        height_arr.append(result.params['p1'].value)

        # create data matrix to fit Dunham coefficients
        hlp = [vg, Jg, ve, Je, result.params['p2'].value, line_info[k]['iso']]

        data.append(hlp)

    else:
        if Jg == 1 and Je == 1:

            (result_multi, xf, yf) = fit_multi_lines(x, y)
            
            plt.plot(xf, yf, '-')
            plt.plot(x, y, '.-')

            plt.text(np.mean(xf), np.max(yf), str(k))
        
            # adding it manually
            q_result = [[vg, Jg, ve, Je, result_multi.params['p2'].value, 35],\
                        [vg, Jg, ve, Je, result_multi.params['p2'].value + result_multi.params['isotope_shift'].value, 37]]
            


#######################################################################
# Fit Dunham coefficients
#######################################################################

Ug = [
[0.0, 3.5429054797521915, -0.00002, 1.0e-7],
[1894.8535149119884, 1],
[1.8535149119884],
]

Ue = [
[38255.25858386326 + wavemeter_offset/100.0/c, 3.716005571623164, -0.00013660808715425417],
[1736.3491645364156, 0.0056375518206914456],
]



(result, initial_guesses) = do_fit(data, Ug, Ue)

print()
print_params(result.params)

(Ug, Ue) = make_dunham_from_params(result.params)



data2losalamos(data, add_q = True, q_lines = q_result, save_filename = 'data_00.csv')


make_report(data, Ug, Ue, vmax = 1, Jmax = 1, save_filename = 'lines_dunham.txt')


####################################
## Predict new lines
####################################
#
#predict_lines(Ug, Ue)
#
#predict_lines(Ug, Ue, vg = 0, ve = 1)

###################################
# Plot complete spectrum
###################################

cnt_freq = 1146.330e12

plt.figure()
for k in range(len(x_arr)):

    x = x_arr[k]
    y = y_arr[k]

    if line_type[k] == 'Q':
        color = 'k'
    elif line_type[k] == 'R':
        color = 'r'
    elif line_type[k] == 'P':
        color = 'b'

    plt.plot((x - cnt_freq)/1e9, 2*60*y, '.-', color = color)

   

(nus, y35, y37) = get_spectrum(Ug, Ue, vg = 0, ve = 0, Jmax = 10)

nus_plot = (nus - cnt_freq)/1e9

my_label = ''

a = 1

plt.plot(nus_plot, a*y35[0], 'r-', label = my_label + '-P')
plt.plot(nus_plot, a*y35[1], 'k-', label = my_label + '-Q')
plt.plot(nus_plot, a*y35[2], 'b-', label = my_label + '-R')

plt.plot(nus_plot, a*y37[0], 'r--')
plt.plot(nus_plot, a*y37[1], 'k--')
plt.plot(nus_plot, a*y37[2], 'b--')

plt.xlim([np.min(nus_plot), np.max(nus_plot)])

plt.ylabel("Signal (a.u.)")

plt.xlabel("Frequency (GHz) + {0:3.6f} THz".format(cnt_freq/1e12))
plt.yticks([])

plt.legend()

plt.tight_layout() 



plt.grid()

plt.show()


