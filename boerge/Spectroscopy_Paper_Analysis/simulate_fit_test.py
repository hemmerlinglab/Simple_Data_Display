import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import sys

sys.path.append("../Analysis_Scripts/")
from constants import *
from data_functions import *
from aux_spectrum_functions import get_spectrum, get_gaussians

from dunham_fit_functions import fit_dunham_coefficients, get_U_init


###############################################

def plot_prediction(x, y, fac = 1.0, offset = 0.0, style = '-', cnt_freq = 0.0, color = None):

    if color is None:
        c1 = 'b'
        c2 = 'k'
        c3 = 'r'
    else:
        c1 = color[0]
        c2 = color[1]
        c3 = color[2]
    
    plt.plot((x - cnt_freq)/1e9, fac*y['P'] + offset, style, color = c1)
    plt.plot((x - cnt_freq)/1e9, fac*y['Q'] + offset, style, color = c2)
    plt.plot((x - cnt_freq)/1e9, fac*y['R'] + offset, style, color = c3)


def plot_simulated_spectrum(line_data, cnt_freq = 0.0):

    df = 100e9
    x = np.linspace(cnt_freq - df, cnt_freq + df, 5000)
    
    y = np.zeros(len(x))

    cnts = list(map(lambda x : x[4], line_data))

    y = get_gaussians(x, cnts, np.ones(len(line_data)), T = 4)

    #plt.plot((x - cnt_freq)/1e9, y)

    for k in range(len(cnts)):
        plt.axvline((cnts[k] - cnt_freq)/1e9, ls = '--', color = 'r')

    plt.xlim( np.min(x - cnt_freq)/1e9, np.max(x - cnt_freq)/1e9 )

    #plt.show()

    return


#########################################################################################################

cnt_freq_00 = 1146.330000e12
cnt_freq_11 = 1145.330000e12 - 100e9


## load fitted data
fit_all = pickle.load( open( 'simulated_data.pickle', "rb" ) )





(Ug_init, Ue_init) = get_U_init()


Ue_init = [
[ 38254.32541337087,3.7363766737869124 ],
[ 1752.0472289516722,-0.1584818467729563 ],
[ -57.895146873435664 ],
[0.1]
]





(Ug, Ue) = fit_dunham_coefficients(fit_all['line_data'], Ug_init, Ue_init, vary_groundstate = False)

comparison = make_report(fit_all['line_data'], Ug, Ue)#, latex_dunham_file = 'latex_dunham.tex', latex_prediction_file = 'latex_prediction.tex')


print_dunham(Ue)

##################################
# Plot errors
##################################

plot_errors(comparison, my_max = -1)

###################################
## Plot spectrum and theory
###################################

(nus_00, y35_00, y37_00) = get_spectrum(Ug, Ue, vg = 0, ve = 0, Jmax = 15, df = 120e9, T = 10)

(nus_11, y35_11, y37_11) = get_spectrum(Ug, Ue, vg = 1, ve = 1, Jmax = 10, df = 120e9, T = 10)


mov_avg = True
mov_n = 5

plt.figure()

plot_simulated_spectrum(fit_all['line_data'], cnt_freq = cnt_freq_00)

plot_prediction(nus_00, y35_00, fac = 2, offset = 0.0, cnt_freq = cnt_freq_00, color = ['k'] * 3)
plot_prediction(nus_00, y37_00, fac = 4, offset = 0.0, style = '--', cnt_freq = cnt_freq_00, color = ['k'] * 3)


plt.figure()

plot_simulated_spectrum(fit_all['line_data'], cnt_freq = cnt_freq_11)

plot_prediction(nus_11, y35_11, fac = 2, offset = 0.0, cnt_freq = cnt_freq_11, color = ['k'] * 3)
plot_prediction(nus_11, y37_11, fac = 4, offset = 0.0, style = '--', cnt_freq = cnt_freq_11, color = ['k'] * 3)



plt.show()



