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


def plot_Ue_error(M, avg_err):

    plt.figure(figsize = (10,8))

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

def plot_Ue_vs_avg_error(M, avg_err):

    plt.figure(figsize = (10,8))

    err = []

    for k in range(len(M[0])):
        for l in range(len(M[0][k])):

            hlp = list(map( lambda x : x[k][l], M ))

            err.append([k, l, hlp])

    no_rows = 3
    no_cols = 2

    for n in range(len(err)):
       
        k = err[n][0]
        l = err[n][1]

        plt.subplot(no_rows, no_cols, k * no_cols + l + 1)

        mean_val = np.mean(err[n][2])
        err_dist = mean_val - err[n][2]

        my_max = np.max(np.abs(err_dist))
        
        scaling = 1

        err_dist /= scaling

        plt.plot( err_dist, avg_err, 'o' )

        plt.xlim(-my_max, my_max)
    
        plt.xlabel('U[{0}][{1}] + {2:.4f} (1/cm)'.format(k, l, mean_val))

        #print(avg_err, err_dist)

    plt.tight_layout()

    return


############################################################################


# load error simulation
arr = pickle.load( open( 'error_simulation.pickle', "rb" ) )


plot_Ue_error(arr['Ue_results'], arr['avg_err_arr'])

plot_Ue_error(arr['De_results'], arr['avg_err_arr'])

plot_Ue_vs_avg_error(arr['Ue_results'], arr['avg_err_arr'])

plot_Ue_vs_avg_error(arr['De_results'], arr['avg_err_arr'])

plt.show()



