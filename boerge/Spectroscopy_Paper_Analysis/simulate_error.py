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


def plot_Ue_error(M, avg_err):

    plt.figure(figsize = (10,6))

    err = []

    for k in range(len(M[0])):
        for l in range(len(M[0][k])):

            hlp = list(map( lambda x : x[k][l], M ))

            err.append([k, l, hlp])

    no_rows = 3
    no_cols = 3

    for n in range(len(err)):
       
        k = err[n][0]
        l = err[n][1]

        plt.subplot(no_rows, no_cols, k * no_cols + l + 1)

        mean_val = np.mean(err[n][2])
        err_dist = mean_val - err[n][2]

        my_max = np.max(np.abs(err_dist))
        
        scaling = 1

        err_dist /= scaling

        plt.hist( err_dist, bins = np.int(len(M)/3.0) )


        plt.xlim(-my_max, my_max)
    
        plt.xlabel('U[{0}][{1}] + {2:.4f} (1/cm)'.format(k, l, mean_val))

    #plt.tight_layout()

    plt.subplot(no_rows, no_cols, no_rows * no_cols)
    plt.hist( avg_err, bins = np.int(len(M)/3.0), color = 'red' )

    plt.xlabel('Avg. Line Error (MHz)')

    plt.tight_layout()

    return



#####################################################################################

# load fitted data
fit_all = pickle.load( open( 'fitted_data.pickle', "rb" ) )

line_data = fit_all['line_data']

#print(fit_all)
#
#asd

#shifted_data = [shifted_data[i] for i in [0,1,2,3,5,6]]#,7,8,9,10,11,12]]


(Ug_init, Ue_init, Dg_init, De_init) = get_U_init()

#Ue_init = [
#    [38255.2403528187, 3.735733648830435], #2.3895804164440537e-06],
#    [1735.1086783255587, -0.15628385798669306],
#    [37.95514890891632],
#    [-1.21762341993954],
#    ]


# Do a Monte Carlo simulation

Ue_results = []
De_results = []
avg_err_arr = []

line_error = 2.5e6

for k in range(500):
    
    if k % 50 == 0:
        print(k)

    # add random shift to data
    shifted_data = add_shifts(line_data.copy(), line_error = line_error)

    (Ug, Ue, Dg, De) = fit_dunham_coefficients(shifted_data, Ug_init, Ue_init, Dg_init, De_init, vary_groundstate = False)

    all_result, avg_err = compare_exp_theory(shifted_data, Ug, Ue, Dg, De)
    
    Ue_results.append(Ue)
    De_results.append(De)
    avg_err_arr.append(avg_err)


plot_Ue_error(Ue_results, avg_err_arr)

plot_Ue_error(De_results, avg_err_arr)


# save simulation results

with open('error_simulation.pickle', 'wb') as f:
        pickle.dump([Ue_results, De_results, avg_err_arr, line_error], f)


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

