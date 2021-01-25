import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import sys

sys.path.append("../Analysis_Scripts/")
from constants import *
from data_functions import *
from energy_functions import *
from aux_spectrum_functions import get_spectrum

from dunham_fit_functions import fit_dunham_coefficients, get_U_init


def calc_lines(Yg, Ye, vg = 0, ve = 0, Jmax = 10, isotope = 35):
    
    line_data = []
    for Jg in range(Jmax):
        for Je in range(Jmax):

            freq = get_energy(Yg, Ye, vg, Jg, ve, Je) 
            
            hlp = [vg, Jg, ve, Je, freq, isotope]
       
            add_data = False
            if (Jg == 1) and (Je == 1):
                add_data = True
            elif (Je == Jg + 1) and (Je >= 1):
                add_data = True
            elif (Je == Jg - 1) and (Je >= 1):
                add_data = True

            if add_data:
                line_data.append(hlp)

    return line_data


cnt_freq_00 = 1146.330000e12
cnt_freq_11 = 1145.330000e12 - 100e9





line_data = []

(Ug_init, Ue_init) = get_U_init()


Ue_init = [
[ 38254.32541337087,3.7363766737869124 ],
[ 1752.0472289516722,-0.1584818467729563 ],
[ -57.895146873435664 ],
]

Jmax = 5

(Yg35, Ye35, Yg37, Ye37) = get_scaled_dunham(Ug_init, Ue_init)

line_data = calc_lines(Yg35, Ye35, Jmax = Jmax, isotope = 35)
line_data.extend(calc_lines(Yg37, Ye37, Jmax = Jmax, isotope = 37))

line_data.extend(calc_lines(Yg35, Ye35, vg = 1, ve = 1, Jmax = Jmax, isotope = 35))

hlp = calc_lines(Yg37, Ye37, vg = 1, ve = 1, Jmax = 2, isotope = 37)

print(hlp)
line_data.append(hlp[1])





# save fit_all dictionary to simulate error
with open('simulated_data.pickle', 'wb') as f:
    pickle.dump({'line_data' : line_data}, f)





