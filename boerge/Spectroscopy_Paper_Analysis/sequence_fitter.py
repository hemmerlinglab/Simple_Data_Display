import sys

sys.path.append("../Analysis_Scripts/")
from constants import *

from dunham_fit_functions import fit_dunham_coefficients, get_U_init
from data_functions import compare_exp_theory



######################################################################

def do_sequence_fit(data):

    my_lines = data[0]['line_data']


    (Ug_init, Ue_init) = get_U_init()


    Ue_init = [
    [ 38254.32541337087,3.7363766737869124 ],
    [ 1752.0472289516722,-0.1584818467729563 ],
    [ -57.895146873435664 ],
    [0.1]
    ]

    Ue_init = [
    [ 38254.32541337087, 3.7363766737869124],# 0.1],#, 0.0001],
    [ 1752.0472289516722, 0.1],
    [ 0.1 ]
    ]

    my_lines = data[0]['line_data'][:]
    
    my_lines.extend(data[1]['line_data'][:])

    #for k in range(len(my_lines)):
    #    print(my_lines[k])

    (Ug, Ue) = fit_dunham_coefficients(my_lines, Ug_init, Ue_init, vary_groundstate = False)

    (all_result, avg_err) = compare_exp_theory(my_lines, Ug, Ue)

    #print(my_lines)
    
    print(avg_err)

    print(Ue[0])
    print(Ue[1])
    #print(Ue[2])

    return


