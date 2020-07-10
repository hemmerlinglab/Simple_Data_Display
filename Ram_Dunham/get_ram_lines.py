from fit_dunham import *
import numpy as np
from helper import get_params_energy, scale_coeff, print_params, make_report, print_matrix, reduce_dunham


c = 299792458

def get_specific_line(params, v1, J1, v2, J2):

    eng = get_params_energy(params, v1, J1, v2, J2)

    print("({0}, {1}) -> ({2}, {3}) = {4:8.2f} = {5:12.6f} THz = {6:12.6f} THz (IR)".format(
        v1, J1, v2, J2, eng, eng * 100.0 * c/1e12, eng * 100.0 * c / 3e12))
    print()

def invert_Y(M):

    Mnew = [ [M[0][0], M[1][0]] ,
            [M[0][1], M[1][1]],
            [M[0][2]],
            [M[0][3]],
            [M[0][4]]]

    massAl = 26.98153841
    massCl_35 = 34.96885269
    massCl_37 = 36.96590258

    U = scale_coeff(Mnew, massAl, massCl_35, scale = False)

    return U
   

# d = [ [gs_v, gs_J, ex_v, ex_J, freq] ]

d = [ 
        [0,0,0,0, 38237.00],
        [1,0,1,0, 38201.72],
        [2,0,2,0, 38157.21],
        [3,0,3,0, 38105.26],
        [4,0,4,0, 38043.44],
        [4,0,5,0, 38433.85],
        [5,0,5,0, 37972.12], # typo in paper
        [6,0,6,0, 37888.98],
        [7,0,6,0, 37435.88], 
        [9,0,7,0, 36898.34],
        [10,0,8,0, 36794.91],
        [15,0,10,0, 35248.74]
        ]

d = [ 
        [0,1,0,1, 38237.00],
        [1,1,1,1, 38201.72],
        [2,1,2,1, 38157.21],
        [3,1,3,1, 38105.26],
        [4,1,4,1, 38043.44],
        [4,1,5,1, 38433.85],
        [5,1,5,1, 37972.12], # typo in paper
        [6,1,6,1, 37888.98],
        [7,1,6,1, 37435.88], 
        [9,1,7,1, 36898.34],
        [10,1,8,1, 36794.91],
        [15,1,10,1, 35248.74]
        ]


# rotational states

d.append([4,40,4,40,38036.01])
d.append([4,48,4,48,38032.55])
d.append([6,30,6,30,37880.80])
d.append([9,25,7,25,36892.90])

d.append([3,42,3,41,38081.28])

d.append([3,51,3,51,38097.80])
d.append([3,62,3,62,38093.83])

d.append([4,40,4,40,38036.01])

d.append([4,63,5,63,38399.35])


# initial guesses

massAl = 26.98153841
massCl_35 = 34.96885269
massCl_37 = 36.96590258

mu = massAl * massCl_35/(massAl + massCl_35)

Ug = [
[0.0, 3.5429054797521915],#, -0.00002],
[1894.8535149119884, -0.0956137436100693],
[-40.15799363940076],
[2.0460078460996782],
[-0.06902618938764313],
]

Ue = [
[38252.04903669513, 3.6943053466270657],#, -0.00002],
[1807.2349511600205, -0.22389431070500337],
[-102.1527235080945],
[2.9407617704423545],
[-2.761287211361239],
]



(result, initial_guesses) = do_fit(d, [Ug, Ue])

#print_params(initial_guesses)
#print_params(result)


print("(( v, J) -> ( v', J')) = <paper> / <prediction> diff: <difference>")

total_err = 0

for n in range(len(d)):

    fit = fcn2min(result, [], [d[n]], get_fit = True)[0]

    err = d[n][4] - fit

    total_err += np.abs(err)

    print('(' + "({0:2d},{1:2d}) -> ({2:2d},{3:2d})".format(d[n][0], d[n][1], d[n][2], d[n][3]) + ') = ' + "{0:8.2f} / {1:8.2f}   diff: {2:8.2f} 1/cm".format(d[n][4], fit, err))


print("\nTotal error/no of lines = {0:6.5f}\n".format(total_err/len(d)))




###############################################
# calculate predictions
###############################################

print()
print_params(result)


(Yg, Ye) = make_params_dunham(result)

#print_matrix(Yg)
#print_matrix(Ye)

(Ug, Ue) = reduce_dunham(Yg, Ye)

print_matrix(Ug)
print_matrix(Ue)
make_report(Ug, Ue, vmax = 1, Jmax = 1, save_filename = 'lines_dunham.txt')


#print_matrix(Ug)
#print_matrix(Ue)

