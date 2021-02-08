import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import sys

sys.path.append("../Analysis_Scripts/")
from constants import *
from data_functions import *
from aux_spectrum_functions import get_spectrum

from dunham_fit_functions import fit_dunham_coefficients, get_U_init

import copy

###############################################
# Change global font size

font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)
###############################################


cnt_freq_00 = 1146.330000e12
cnt_freq_11 = 1145.330000e12 - 100e9




def add_shifts(d, line_error = 5.0e6):

    # this function adds random shifts to each measured frequency

    for k in range(len(d)):
        #d[k][4] += (2*np.random.rand() - 1.0)*10e6
        
        d[k][4] += np.random.normal(loc = 0.0, scale = line_error)

    return d


#####################################################################################

# load fitted data
fit_all = pickle.load( open( 'fitted_data.pickle', "rb" ) )

line_data = fit_all['line_data']

#print(fit_all)



(Ug_init, Ue_init, Dg_init, De_init) = get_U_init()


# Do a Monte Carlo simulation

Ue_results = []
De_results = []
avg_err_arr = []

line_error = 10.0e6


for k in range(5000):
    
    if k % 50 == 0:
        print(k)

    # add random shift to data

    # a deep copy is required here to REALLY copy the lsit
    hlp_line_data = copy.deepcopy(line_data)

    shifted_data = add_shifts(hlp_line_data, line_error = line_error)
    
    (Ug, Ue, Dg, De) = fit_dunham_coefficients(shifted_data, Ug_init, Ue_init, Dg_init, De_init, vary_groundstate = False)

    all_result, avg_err = compare_exp_theory(shifted_data, Ug, Ue, Dg, De)
    
    Ue_results.append(Ue)
    De_results.append(De)
    avg_err_arr.append(avg_err)


#plot_Ue_error(Ue_results, avg_err_arr)

#plot_Ue_error(De_results, avg_err_arr)


# save simulation results

with open('error_simulation.pickle', 'wb') as f:
    pickle.dump({'Ue' : Ue, 'Ue_results' : Ue_results, 'De' : De, 'De_results' : De_results, 'avg_err_arr' : avg_err_arr, 'line_error' : line_error}, f)


plt.show()




#comparison = make_report(line_data, Ug, Ue, Dg, De, latex_dunham_file = 'latex_dunham.tex', latex_prediction_file = 'latex_prediction.tex')
#
#plot_errors(comparison)
#
#
#(Yg35, Ye35, Yg37, Ye37) = get_scaled_dunham(Ug, Ue, Dg, De)

#print_dunham(Yg35)
#print_dunham(Yg37)
#print_dunham(Ue)

#print(Ug)
#print(Ue)
#print(avg_err)


plt.show()

