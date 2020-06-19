import numpy as np
import lmfit
from lmfit import Minimizer, Parameters, report_fit
from prettytable import PrettyTable



# global constants

c = 299792458
kB = 1.3806482e-23
h_planck = 6.63607e-34
amu = 1.66e-27

massAl = 26.98153841
massCl_35 = 34.96885269
massCl_37 = 36.96590258


#############################################################

def print_params(params):

    print()
    for k in sorted(params.keys()):

        par = params[k]

        print("{0} : {2:12.6f} <= {1:12.6f} <= {3:12.6f}".format(k, par.value, par.min, par.max))


def print_matrix(U, txt = None):

    if not txt is None:
        print(txt + ' = [')
    else:
        print('[')
    for k in range(len(U)):
        print(U[k], end = '')
        #for l in range(len(U[k])):
        #    print(U[k][l])
        if k < len(U):
            print(',')

    print(']\n')


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

    return (Ug, Ue)



def save_dunham(Ug, Ue, ext = ''):

    # save coefficients
    f = open('Ug_' + 'ram' + '.txt', 'w')
    for k in range(len(Ug)):
        f.write(",".join(map(str, Ug[k])))
        f.write("\n")
    f.close()

    f = open('Ue_' + 'ram' + '.txt', 'w')
    for k in range(len(Ue)):
        f.write(",".join(map(str, Ue[k])))
        f.write("\n")
    f.close()


def load_dunham(Ug_file, Ue_file):

    # read in Dunham matrices
    f = open(Ug_file)
    lines = f.readlines()
    Ug = []
    for line in lines:
        Ug.append( [float(i) for i in line.strip().split(',')] )
    f.close()
    
    f = open(Ue_file)
    lines = f.readlines()
    Ue = []
    for line in lines:
        Ue.append( [float(i) for i in line.strip().split(',')] )
    f.close()

    return (Ug, Ue)



def scale_dunham_matrix(M, mass1, mass2, scale = True):

    # scale or unscale Dunham coefficients with the reduced mass
    mu = (mass1 * mass2)/(mass1 + mass2)

    if scale == False:
        
        # return U_kl = 1/mu^(-k/2 - l) * M_kl
        U = []
        for k in range(len(M)):
            hlp = []
            for l in range(len(M[k])):
                hlp.append(1.0/mu**(-k/2.0 - l) * M[k][l])
            U.append(hlp)
 
        return U
    else:

        # return Y_kl = mu^(-k/2 - l) * M_kl
        Y = []
        for k in range(len(M)):
            hlp = []
            for l in range(len(M[k])):
                hlp.append(mu**(-k/2.0 - l) * M[k][l])
            Y.append(hlp)
        
        return Y



def get_dunham(Ug, Ue):

    
    # ground state
    Yg35 = scale_dunham_matrix(Ug, massAl, massCl_35, scale = True)
    Yg37 = scale_dunham_matrix(Ug, massAl, massCl_37, scale = True)

    # excited state
    Ye35 = scale_dunham_matrix(Ue, massAl, massCl_35, scale = True)
    Ye37 = scale_dunham_matrix(Ue, massAl, massCl_37, scale = True)

    return (Yg35, Ye35, Yg37, Ye37)


def reduce_dunham35(Yg, Ye, isotope = 35):

    if isotope == 35:
        mass2 = massCl_35
    elif isotope == 37:
        mass2 = massCl_37
    else:
        print('Error')
        asd

    Ug = scale_dunham_matrix(Yg, massAl, massCl_35, scale = False)
    Ue = scale_dunham_matrix(Ye, massAl, massCl_35, scale = False)

    return (Ug, Ue)

def get_line_type(Jg, Je):
    
    if Je == Jg:
        my_type = 'Q'
    elif Je == Jg - 1:
        my_type = 'P'
    elif Je == Jg + 1:
        my_type = 'R'

    return my_type

def make_report(data, Ug, Ue, vmax = 1, Jmax = 1, save_filename = 'report.txt'):

    
    (Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)

    print()
    print("*"*80)
    print("Report")
    print("*"*80)

    # resulting dunham coefficients
    print("\n" + "*"*40)
    print("X-state")
    print("*"*40 + "\n")

    print_matrix(Ug, txt = 'Ug')
    print_matrix(Yg35, txt = 'Yg-35')
    print_matrix(Yg37, txt = 'Yg-37')

    print()

    # resulting dunham coefficients
    print("\n" + "*"*40)
    print("A-state")
    print("*"*40 + "\n")

    print_matrix(Ue, txt = 'Ue')
    print_matrix(Ye35, txt = 'Ye-35')
    print_matrix(Ye37, txt = 'Ye-37')

    print()
    save_dunham(Ug, Ue, '_fit')
   

    # table to compare measured and calculated lines
    
    t = PrettyTable(['vg', 'Jg', 've', 'Je', 'Exp. Freq. (THz)', 'Cal. Freq. (THz)', 'Diff. (MHz)', 'AlCl', 'Type'])

    t.float_format['Exp. Freq. (THz)'] = ".6"
    t.float_format['Cal. Freq. (THz)'] = ".6"
    t.float_format['Diff. (MHz)'] = "8.2"

    total_err = 0.0

    for k in range(len(data)):
        
        # d = [ [gs_v, gs_J, ex_v, ex_J, freq, 35/37] ]

        vg = data[k][0]
        Jg = data[k][1]
        ve = data[k][2]
        Je = data[k][3]
        meas_freq = data[k][4]
        isotope = data[k][5]

        
        if isotope == 35:
            eng = get_energy(Yg35, Ye35, vg, Jg, ve, Je)
        elif isotope == 37:
            eng = get_energy(Yg37, Ye37, vg, Jg, ve, Je)
        else:
            print('Error')
            asd

        my_type = get_line_type(Jg, Je) 
        
        t.add_row([vg, Jg, ve, Je, meas_freq/1e12, eng/1e12, (eng - meas_freq)/1e6, isotope, my_type])

        total_err += np.abs((eng - meas_freq)/1e6)


    t.add_row(['-']*6 + ['--------'] + ['-']*2)
    t.add_row(['','','','','','','avg:' + "{0:2.2f}".format(total_err/len(data)),'',''])
    
    print(t)
    


def predict_lines(Ug, Ue, save_filename = 'line_predictions', vg = 0, ve = 0, Jmax = 10):

    (Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)
    
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
 
        if my_type is 'P':
            color = '34'
        elif my_type is 'Q':
            color = '32'
        elif my_type is 'R':
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

    print(t)


    f = open(save_filename, 'w')
    np.savetxt(f, d, delimiter = ',')
    f.close()

    # print order-reduced coefficients

    # (we  + wexe  * (v+1/2) + weye * (v+1/2)^2) * (v+1/2)
    # (we' + wexe' * (v+1/2)) * (v+1/2)
    # (we'') * (v+1/2)



    return 





##########################################################################


def scale_dunham(U, k, l, isotope):

    mass1 = massAl # aluminum
    
    if isotope == 35:
        mass2 = massCl_35
    elif isotope == 37:
        mass2 = massCl_37
    
    mu = (mass1 * mass2)/(mass1 + mass2)

    return (mu**(-k/2.0 - l) * U)


def calc_energy(val, v, J, k, l):
    return 100 * c * (val * (v + 0.5)**k * ( J * (J+1.0) )**l)


def energy(Y, v, J):

    e = 0.0
    for k in range(len(Y)):
        for l in range(len(Y[k])):

            e += calc_energy(Y[k][l], v, J, k, l)

    return e


def get_energy(Yg, Ye, vg, Jg, ve, Je):
    
    return energy(Ye, ve, Je) - energy(Yg, vg, Jg)



def get_params_energy(params, vg, Jg, ve, Je, isotope):

    e = 0.0
    for key in params.keys():

        # Ue35
        k = int(key[2])
        l = int(key[3])

        state = key[0:2]
        
        # scale mass-reduced coefficients to the correct isotope
        val = scale_dunham(params[key].value, k, l, isotope)
        
        if state == 'Ue':
            e += calc_energy(val, ve, Je, k, l)

        if state == 'Ug':
            e -= calc_energy(val, vg, Jg, k, l)

    return e



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


def get_params(Ug, Ue):

        # make parameters
        params = Parameters()

        # ground state
        for k in range(len(Ug)):
            for l in range(len(Ug[k])):

                val = Ug[k][l]

                params.add('Ug' + str(k) + str(l), value = val, vary = True)
               
        # excited state
        for k in range(len(Ue)):
            for l in range(len(Ue[k])):

                val = Ue[k][l]
                
                params.add('Ue' + str(k) + str(l), value = val, vary = True)

        # electronic constants
        params['Ug00'].vary = False

        # Use ground state vibrational energy levels from Bernath
        params['Ug10'].vary = False
        params['Ug10'].value = 1880.204
        
        params['Ug20'].vary = False
        params['Ug20'].value = -32.012

        params['Ug11'].vary = False
        params['Ug11'].value = -9.575174e-2

        params['Ug01'].vary = False
        params['Ug01'].value = 3.715174

        params['Ug02'].vary = False
        params['Ug02'].value = -5.802349e-5

        params['Ug03'].vary = False
        params['Ug03'].value = -1.964e-10


        return params



def do_fit(d, Ug, Ue):
 
        params = get_params(Ug, Ue)

        # do fit, here with leastsq model
        minner = Minimizer(fcn2min_dunham, params, fcn_args=([], d))
        result = minner.minimize()
        

        # Store the Confidence data from the fit
        #con_report = lmfit.fit_report(result.params)
        
        #report_fit(result)

        return (result, params)




