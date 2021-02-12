import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit
import sys

sys.path.append("../Analysis_Scripts/")
from constants import *
from energy_functions import *

###################################################################

def get_Ug_Bernath():

    ## unconstrained values from http://dx.doi.org/10.1063/1.465611
    #Ug = [
    #    [0.0, 3.71517423, -5.802349e-05, -1.964e-10, 6.208e-14],
    #    [1880.20433, -9.5757548e-2, 4.0108e-7],
    #    [-32.01271, 1.088211e-3, 2.5939e-8],
    #    [3.95499e-1],
    #    [-4.861e-3]
    #]

    # constrained 
    Ug = [
        [0.0, 3.71517408, -5.80214265e-05, -1.57686816e-10, -0.35748725e-14],
        [1880.20216, -9.5756544e-2, 3.87217711e-7],
        [-32.01210, 1.087167e-3, 2.84009761e-8],
        [3.95186e-1],
        [-4.802e-3]
    ]

    # constrained 
    Dg = [
        [0.0, -1.4432, 0.0, 0.0, 0.0],
        [-1.2238, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0],
        [0.0]
    ]

    return (Ug, Dg)

###################################################################

def get_U_init(set_bob_zero = False):

    (Ug, Dg) = get_Ug_Bernath()

    Ue = [
    [ 38257.17501435353,3.7373561499782517 ],
    [ 1722.3811711416483,-0.16164221289477823 ],
    [ 0.1 ]
    ]

    De =  [
    [ 0.0, 0.0 ],
    [ -1.0, 0.0 ],
    [ 0.0 ]
    ]

    ###################

    #Ue = [
    #[ 38257.17501435353,3.7373561499782517 ],
    #[ 1722.3811711416483,-0.16164221289477823 ],
    #[ 0.1 ]
    #]

    Ue = [
    [ 38253.16953987888,3.7370100040152563 ],
    [ 1766.3041492863101,-0.1609525199616741 ],
    [ -85.77292247914849 ],
    ]

    De =  [
    [ -0.36, 0.0 ],
    [ 0.0, 0.0 ],
    [ 0.0 ]
    ] 

    if set_bob_zero:
        De = [
          [ 0.0, 0.0 ],
          [ 0.0, 0.0 ],
          [ 0.0 ]
          ]

    return (Ug, Ue, Dg, De)



###################################################################

def make_index_list(indeces):

    full_indeces = []

    k_indeces = np.unique(list(map(lambda x : x[0], indeces)))

    for n1 in range(len(k_indeces)):

        full_indeces.append([k_indeces[n1], []])
        
        for n2 in range(len(indeces)):
            # find all l indeces for each k index

            if k_indeces[n1] == indeces[n2][0]:
                full_indeces[n1][1].append(indeces[n2][1])

    return full_indeces

###################################################################

def populate_Y(par, state, ind):

    Y = []
    for n1 in range(len(ind)):
            hlp = []
        
            k = ind[n1][0]
            for n2 in range(len(ind[n1][1])):
                l = ind[n1][1][n2]

                val = par[state + str(k) + str(l)].value
                hlp.append(val)
            
            Y.append(hlp)
            
    return Y

###################################################################

def make_dunham_from_params(par):

    indeces_g = []
    indeces_e = []
    for key in sorted(par.keys()):
        
        k = int(key[2])
        l = int(key[3])

        state = key[0:2]
        if state == 'Ue':
            indeces_e.append([k, l])
        if state == 'Ug':
            indeces_g.append([k, l])

    # get all vibrational indeces

    #k_indeces_e = np.unique(list(map(lambda x : x[0], indeces_e)))

    full_indeces_g = make_index_list(indeces_g)
    full_indeces_e = make_index_list(indeces_e)

    Ug = populate_Y(par, 'Ug', full_indeces_g)
    Ue = populate_Y(par, 'Ue', full_indeces_e)

    Dg = populate_Y(par, 'Dg', full_indeces_g)
    De = populate_Y(par, 'De', full_indeces_e)

    return (Ug, Ue, Dg, De)

###################################################################

def fcn2min_dunham(params, x, data, get_fit = False):

    # data = [[v1, J1, v2, J2, freq], ...]
    model = []
    
    if get_fit == False:

        for d in data:
            model.append( get_params_energy(params, d[0], d[1], d[2], d[3], d[5]) - d[4] )
        
        return np.array(model)

    else:

        for d in data:
            model.append( get_params_energy(params, d[0], d[1], d[2], d[3], d[5]) )

        return np.array(model)

###################################################################

def get_params(Ug, Ue, Dg, De, vary_groundstate):

        # make parameters
        params = Parameters()

        if vary_groundstate:
            
            print("Don't vary groundstate")
            asd
            
            # ground state
            for k in range(len(Ug)):
                for l in range(len(Ug[k])):

                    val = Ug[k][l]

                    params.add('Ug' + str(k) + str(l), value = val, vary = True)
                
                    # adding Born-Oppenheimer correction
                    params.add('Dg' + str(k) + str(l), value = 0.0, vary = False)


        else:

            # use Bernath values and don't vary the values
            Ug, Dg = get_Ug_Bernath()

            for k in range(len(Ug)):
                for l in range(len(Ug[k])):
                    
                    val = Ug[k][l]
                    
                    params.add('Ug' + str(k) + str(l), value = val, vary = False)
    
                    # adding Born-Oppenheimer correction
                    val = Dg[k][l]

                    params.add('Dg' + str(k) + str(l), value = val, vary = False)

        # excited state
        for k in range(len(Ue)):
            for l in range(len(Ue[k])):

                val = Ue[k][l]
                
                params.add('Ue' + str(k) + str(l), value = val, vary = True)

                # Born-Oppenheimer corrections for excited state

                # Should throw error if Ue and De don't have the same dimensions
                
                val = De[k][l]

                # adding Born-Oppenheimer correction

                if not val == 0.0:
                    params.add('De' + str(k) + str(l), value = val, vary = True)
                else:
                    params.add('De' + str(k) + str(l), value = 0.0, vary = False)

        return params

###################################################################

def fit_dunham_coefficients(d, Ug, Ue, Dg, De, vary_groundstate = False):
 
        # get the params from reduced mass matrix
        params = get_params(Ug, Ue, Dg, De, vary_groundstate)

        # do fit, here with leastsq model
        minner = Minimizer(fcn2min_dunham, params, fcn_args=([], d))
        result = minner.minimize()
        
        # Store the Confidence data from the fit
        #con_report = lmfit.fit_report(result.params)
        
        #report_fit(result)

        # transform params back to Dunham
        (Ug, Ue, Dg, De) = make_dunham_from_params(result.params)

        return (Ug, Ue, Dg, De)


