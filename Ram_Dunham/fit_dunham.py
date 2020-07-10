import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit
from helper import get_params_energy, get_dunham, print_params


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

    
def make_params_dunham(par):

    indeces_g = []
    indeces_e = []
    for key in sorted(par.keys()):
        
        k = int(key[2])
        l = int(key[3])

        state = key[0:2]
        if state == 'Ye':
            indeces_e.append([k, l])
        if state == 'Yg':
            indeces_g.append([k, l])

    # get all vibrational indeces

    #k_indeces_e = np.unique(list(map(lambda x : x[0], indeces_e)))

    full_indeces_e = make_index_list(indeces_e)
    full_indeces_g = make_index_list(indeces_g)

    Yg = populate_Y(par, 'Yg', full_indeces_g)
    Ye = populate_Y(par, 'Ye', full_indeces_e)

    return (Yg, Ye)


def fcn2min(params, x, data, get_fit = False):

    a = params['Yg10']
    b = params['Ye00']
    c = params['Ye10']
    
    # data = [[v1, J1, v2, J2, freq], ...]
    model = []
    
    if get_fit == False:

        for d in data:
            model.append( get_params_energy(params, d[0], d[1], d[2], d[3]) - d[4] )
            
            #model.append( (a * d[0]**2 + b * d[0] + c) - d[4] )
       
        #print_params(params)
        #print(model)

        return np.array(model)

    else:

        for d in data:
            #model.append( energy(Ye, d[2], d[3]) - energy(Yg, d[0], d[1]) )
            model.append( get_params_energy(params, d[0], d[1], d[2], d[3]) )

        return np.array(model)



def get_minmax(val):

    pc = 0.05

    if val == 0.0:
        return (-10, 10)
    elif val > 0.0:
        return ((1-pc) * val, (1+pc) * val)
    else:
        return ((1+pc) * val, (1-pc) * val)



def get_params(Yg, Ye):

        # make parameters
        params = Parameters()
        
        for k in range(len(Yg)):
            for l in range(len(Yg[k])):

                val = Yg[k][l]

                (mymin, mymax) = get_minmax(val)

                #params.add('Yg' + str(k) + str(l), value = val, min = mymin, max = mymax, vary = True)
                params.add('Yg' + str(k) + str(l), value = val, vary = True)
               
                val = Ye[k][l]

                (mymin, mymax) = get_minmax(val)

                #params.add('Ye' + str(k) + str(l), value = val, min = mymin, max = mymax, vary = True)
                params.add('Ye' + str(k) + str(l), value = val, vary = True)

        # electronic constants
        params['Yg00'].vary = False

        #print_params(params)

        return params


def do_fit(d, init_guesses, lines = 35):
 
        # init_guesses = [Ug, Ue] = Dunham coefficients        
        (Yg35, Ye35, Yg37, Ye37) = get_dunham(init_guesses[0], init_guesses[1])

        params = get_params(Yg35, Ye35)

        #print_params(params)

        #params['Yg11'].min = -1
        #params['Yg11'].max = 1
        #
        #params['Yg20'].min = -1
        #params['Yg20'].max = 1


        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=([], d))
        result = minner.minimize()
        
        

        # Store the Confidence data from the fit
        #con_report = lmfit.fit_report(result.params)
        
        #report_fit(result)

        #print(params)
        #print(result.params)

        return (result.params, params)




