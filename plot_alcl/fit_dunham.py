import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit


nvib = 5
nrot = 2


def energy(Y, v, J):

    e = 0
    for k in range(nvib):
       l = 0
       e += Y[0][k] * (v + 0.5)**k * ( J * (J + 1.0) )**l
  
    for k in [0,1]:
        l = 1
        e += Y[1][k] * (v + 0.5)**k * ( J * (J + 1.0) )**l

    return e


def make_dunham(params, state = 'Yg'):

    Y = []
    hlp = []
    for k in range(nvib):
        l = 0
        hlp.append(params[state + str(k) + str(l)])
    Y.append(hlp)

    hlp = []
    for k in [0,1]:
        l = 1
        hlp.append(params[state + str(k) + str(l)])

    Y.append(hlp)

    return Y

# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, get_fit = False):
    
    Yg = make_dunham(params, 'Yg')
    Ye = make_dunham(params, 'Ye')

    # data = [[v1, J1, v2, J2, freq], ...]
    model = []


    if get_fit == False:

        for d in data:
            model.append( energy(Ye, d[2], d[3]) - energy(Yg, d[0], d[1])  - d[4] )

        return np.array(model)

    else:

        for d in data:
            model.append( energy(Ye, d[2], d[3]) - energy(Yg, d[0], d[1]) )

        return np.array(model)


def do_fit(d):
        params = Parameters()
 
        for k in range(nvib):
            l = 0
            params.add('Yg' + str(k) + str(l), value=0.0, vary = True) # ground state
            params.add('Ye' + str(k) + str(l), value=0.0, vary = True) # excited state

        params.add('Yg01', value=0.0, vary = True) # ground state
        params.add('Ye01', value=0.0, vary = True) # excited state

        params.add('Yg11', value=0.0, vary = True) # ground state
        params.add('Ye11', value=0.0, vary = True) # excited state

    

        # Y_nvib_nrot
        # Y00 = Te
        # Y10 = we
        # Y20 = -wexe
        # Y30 = weye
        # Y40 = weze
        
        # Y01 = Be
        # Y11 = -alphae
        # Y21 = gammae
        # Y02 = -De
        # Y12 = -betae
        # Y03 = He
        # Y04 = Le

        
        # electronic constants
        params['Yg00'].value =     0.0
        params['Yg00'].vary = False
        params['Ye00'].value = 38255.968

        # vibrational constants
        params['Yg10'].value = 481.4956
        params['Ye10'].value = 449.6235

        params['Yg20'].value = -1.751
        params['Ye20'].value = -6.359

        params['Yg30'].value = 2.5e-7
        params['Ye30'].value = 0.3999

        params['Yg40'].value = -0.00109
        params['Ye40'].value = -0.0420


        ## rotational constants
        params['Yg01'].value = 0.234
        params['Ye01'].value = 0.234

        params['Yg11'].value = -0.001611
        params['Ye11'].value = -0.002519

        # change the min/max of each parameter
        my_range = 1.0
        for k in range(nvib):
            l = 0
            for par in ['Yg' + str(k) + str(l), 'Ye' + str(k) + str(l)]:

                v = params[par].value

                if not (v == 0):
                    params[par].min = (1 - np.sign(v) * my_range) * v
                    params[par].max = (1 + np.sign(v) * my_range) * v
                else:
                    params[par].min = -1.0
                    params[par].max = +1.0

        for k in [0,1]:
            l = 1
            for par in ['Yg' + str(k) + str(l), 'Ye' + str(k) + str(l)]:

                v = params[par].value

                if not (v == 0):
                    params[par].min = (1 - np.sign(v) * my_range) * v
                    params[par].max = (1 + np.sign(v) * my_range) * v
                else:
                    params[par].min = -1.0
                    params[par].max = +1.0



        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=([], d))
        result = minner.minimize()
        
        # Store the Confidence data from the fit
        con_report = lmfit.fit_report(result.params)
        
        return (result.params, params)




