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
###############################################

font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)

###############################################



def predict_lines(save_filename = 'line_predictions.txt', vg = 0, ve = 0, Jmax = 10):

    (Ug, Ue, Dg, De) = load_dunham()

    (Yg35, Ye35, Yg37, Ye37) = get_scaled_dunham(Ug, Ue, Dg, De)
    
    d = []
    for Jg in range(0, Jmax):
        for Je in range(0, Jmax):

            if (np.abs(Je - Jg) <= 1) and not (Je == 0):
               eng35 = get_energy(Yg35, Ye35, vg, Jg, ve, Je)
               eng37 = get_energy(Yg37, Ye37, vg, Jg, ve, Je)

               hlp1 = [vg, Jg, ve, Je, eng35/100.0/c, eng35/1e12, eng35/3e12, 35]
               hlp2 = [vg, Jg, ve, Je, eng37/100.0/c, eng37/1e12, eng37/3e12, 37]

               d.append(hlp1)
               d.append(hlp2)

    d = np.array(d)

    # sort according to frequency

    d = sorted(d, key = lambda d_entry: d_entry[5]) 

    # print out table
    t = PrettyTable(['vg', 'Jg', 've', 'Je', 'Calc. (1/cm)', 'Calc. (THz)', 'Calc. (IR, THz)', 'AlCl', 'Type'])

    t.float_format['Calc. (THz)'] = ".6"
    t.float_format['Calc. (IR, THz)'] = ".6"
    t.float_format['Calc. (1/cm)'] = ".3"
    t.float_format['AlCl'] = ".0"


    for n in d:

        Jg = n[1]
        Je = n[3]

        my_type = get_line_type(Jg, Je)
       
        row = list(n)
        row.append(my_type)
 
        if my_type == 'P':
            color = '34'
        elif my_type == 'Q':
            color = '32'
        elif my_type == 'R':
            color = '31'
        
        row[0] = "{0:1.0f}".format(row[0])
        row[1] = "{0:1.0f}".format(row[1])
        row[2] = "{0:1.0f}".format(row[2])
        row[3] = "{0:1.0f}".format(row[3])
        row[4] = "{0:.3f}".format(row[4])
        row[5] = "{0:3.6f}".format(row[5])
        row[6] = "{0:3.6f}".format(row[6])
        row[7] = "{0:1.0f}".format(row[7])

        for k in range(len(row)):
            row[k] = "\033[1;" + color + 'm' + row[k] + "\033[1;0m"
        

        t.add_row(row)

    print()
    print(t)
    print()

    f = open(save_filename, 'w')
    np.savetxt(f, d, delimiter = ',')
    f.close()

    # print order-reduced coefficients

    # (we  + wexe  * (v+1/2) + weye * (v+1/2)^2) * (v+1/2)
    # (we' + wexe' * (v+1/2)) * (v+1/2)
    # (we'') * (v+1/2)



    return 



#################################################################

def load_dunham(filename = 'dunham_matrices_fit.pickle'):
    
    arr = pickle.load( open( filename, "rb" ) )

    return (arr['Ug'], arr['Ue'], arr['Dg'], arr['De'])


##############################
# Main
##############################



predict_lines(save_filename = 'line_predictions.txt', vg = 0, ve = 0, Jmax = 8)
predict_lines(save_filename = 'line_predictions.txt', vg = 1, ve = 1, Jmax = 8)

predict_lines(save_filename = 'line_predictions.txt', vg = 2, ve = 2, Jmax = 3)
#predict_lines(save_filename = 'line_predictions.txt', vg = 3, ve = 3, Jmax = 4)

predict_lines(save_filename = 'line_predictions.txt', vg = 1, ve = 0, Jmax = 3)

predict_lines(save_filename = 'line_predictions.txt', vg = 0, ve = 1, Jmax = 3)

predict_lines(save_filename = 'line_predictions.txt', vg = 0, ve = 2, Jmax = 2)
predict_lines(save_filename = 'line_predictions.txt', vg = 1, ve = 2, Jmax = 2)
predict_lines(save_filename = 'line_predictions.txt', vg = 1, ve = 3, Jmax = 2)


