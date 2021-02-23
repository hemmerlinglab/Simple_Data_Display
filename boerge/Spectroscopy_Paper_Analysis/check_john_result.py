import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import sys

sys.path.append("../Analysis_Scripts/")
from constants import *
from data_functions import *
from aux_spectrum_functions import get_transition_energies, get_spectrum

from dunham_fit_functions import fit_dunham_coefficients, get_U_init, get_U_John

###############################################
# Change global font size

font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)
###############################################



# load fitted data
fit_all = pickle.load( open( 'fitted_data.pickle', "rb" ) )

line_data = fit_all['line_data']


(Ug, Ue, Dg, De) = get_U_init()
#(Ug, Ue, Dg, De) = get_U_John()



#all_result, avg_err = compare_exp_theory(line_data, Ug, Ue, Dg, De)


comparison = make_report(line_data, Ug, Ue, Dg, De)

plot_errors(comparison)



plt.show()


