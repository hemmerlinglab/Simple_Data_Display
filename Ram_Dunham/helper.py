import numpy as np
from configparser import ConfigParser
import ast
import matplotlib.pyplot as plt

c = 299792458
kB = 1.3806482e-23
h_planck = 6.63607e-34
amu = 1.66e-27

def get_pop(v, J, we, B, T):

    # returns population distribution

    Ewe = we * h_planck
    EB = B * h_planck

    beta = 1/(kB*T)

    Ev = Ewe * (v+0.5)
    
    Erot = lambda J : EB * J * (J+1)

    # rotational distribution
    J_arr = np.arange(0, 200, 1)

    #print(Jarr)

    Nrot = np.sum( (2*J_arr+1) * np.exp(-beta * Erot(J_arr)) )
    
    pop_rot = 1/Nrot * (2*J+1) * np.exp(-beta * Erot(J))

    # vibrational distribution
    Nvib = np.exp(beta * Ewe/2)/(np.exp(beta * Ewe) - 1)

    pop_vib = 1/Nvib * np.exp(-beta * Ev)

    return pop_vib * pop_rot


#def get_pop(J, B, T, N0 = 1):
#def get_pop(v, J, we, B, T, N0 = 1):
#    # returns rotational population distribution
#
#    pop = N0 * (2*J+1) * np.exp(-B * h_planck * J*(J+1)/(kB*T))
#
#    return pop / (np.sqrt(kB*T/(2*h_planck*B)) - 1.0/2.0)


def get_doppler(T, f0):
    # returns Doppler width for 27Al35Cl

    return np.sqrt(8*kB*T*np.log(2)/((35+27)*amu*c**2)) * f0


def simple_gauss(x, x0, A, w):
    return A * np.exp( -(x-x0)**2/(2*w**2) )


def get_line_centers(Yg, Ye, ve, vg, Jmax = 1, T = 50, real_amplitudes = True):

    Jg_arr = np.arange(0, Jmax+1)    
    Je_arr = np.arange(1, Jmax+1)
    
    f_P = []
    f_Q = []
    f_R = []
    
    A_P = []
    A_Q = []
    A_R = []
    
    for Jg in Jg_arr:
        for Je in Je_arr:
            
            # only dipole transitions
            eng = 100*c*(energy(Ye, ve, Je) - energy(Yg, vg, Jg))
            
            # apply population of ground states
            if real_amplitudes:
                A = get_pop(vg, Jg, 100*c*Yg[1][0], 100*c*Yg[0][1], T)
            else:
                A = 1.0

            if Je - Jg == -1: 
                f_P.append(eng)
                A_P.append(A)
    
            if Je - Jg ==  0: 
                f_Q.append(eng)
                A_Q.append(A)
    
            if Je - Jg == +1: 
                f_R.append(eng)
                A_R.append(A)
    
    f_P = np.array(f_P)
    f_Q = np.array(f_Q)
    f_R = np.array(f_R)
    
    A_P = np.array(A_P)
    A_Q = np.array(A_Q)
    A_R = np.array(A_R)

    return ([f_P, f_Q, f_R], [A_P, A_Q, A_R])


def scale_coeff(M, mass1, mass2, scale = True):

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

def get_params_energy(params, vg, Jg, ve, Je):

    e = 0.0
    for key in params.keys():

        # Ye35
        k = int(key[2])
        l = int(key[3])

        state = key[0:2]
        
        if state == 'Ye':
            val = params[key].value
            e += val * (ve + 0.5)**k * ( Je * (Je + 1.0) )**l

        if state == 'Yg':
            val = params[key].value
            e -= val * (vg + 0.5)**k * ( Jg * (Jg + 1.0) )**l

    return e


def energy(Y, v, J):

    e = 0.0
    for k in range(len(Y)):
        for l in range(len(Y[k])):

            e += Y[k][l] * (v + 0.5)**k * ( J * (J + 1.0) )**l

    return e



def get_energy(Yg, Ye, vg, Jg, ve, Je):
    
    return energy(Ye, ve, Je) - energy(Yg, vg, Jg)



def get_reduced_dunham():

    # ground state from Bernath
    Ug = [
            #[0.0,3.59517408,-5.802142e-5],
            [0.0,3.711517408,-5.802142e-5],
            [1880.202,-9.575654e-2],
            [-32.012],
            [3.95186e-1],
            [-4.802e-3]
            ]
    
    # excited state
    Ue = [
           #[38251.101190451634 - (39.6e9+15.7e9-55.3e9)/100/c, 3.6400046656060443, 0.00], 
           [38251.101190451634 - (39.6e9+15.7e9-55.3e9)/100/c, 3.695, 0.00], 
           #[38237.483145, 3.64, 0.00], 
           [1784.3665982617565, 0.0*-0.2344361133562231],
           [-114.5238709407005], 
           [17.527497484431063], 
           [-6.316749292366754]
         ]

    Ug = [[0, 3.697689639, -0.00004110220051, 0, 0, 0], [1880.20433, \
0.04887253993, 0, 0, 0, 0], [-32.01271, 0, 0, 0, 0, 0], [0.0395499, \
0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

    Ue = [[38255.63489, 3.718219212, 0.00273806453, 0, 0, 0], [1729.598167, \
-0.06421374131, 0, 0, 0, 0], [60.74391688, 0, 0, 0, 0, 0], \
[-180.1417929, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, \
0]]

    return (Ug, Ue)

def get_dunham(Ug, Ue):

    massAl = 26.98153841
    massCl_35 = 34.96885269
    massCl_37 = 36.96590258

    # ground state
    Yg35 = scale_coeff(Ug, massAl, massCl_35, scale = True)
    Yg37 = scale_coeff(Ug, massAl, massCl_37, scale = True)

    # excited state
    Ye35 = scale_coeff(Ue, massAl, massCl_35, scale = True)
    Ye37 = scale_coeff(Ue, massAl, massCl_37, scale = True)

    return (Yg35, Ye35, Yg37, Ye37)


def reduce_dunham(Yg, Ye):

    massAl = 26.98153841
    massCl_35 = 34.96885269
    massCl_37 = 36.96590258

    Ug = scale_coeff(Yg, massAl, massCl_35, scale = False)
    Ue = scale_coeff(Ye, massAl, massCl_35, scale = False)

    return (Ug, Ue)


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
   


# Cl-35 and Cl-37 lines
# summed coefficients for less orders
# plots


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




def make_report(Ug, Ue, vmax = 1, Jmax = 1, save_filename = 'report.txt'):

    
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
    save_dunham(Ug, Ue, '_ram')
   

    vmax += 1
    Jmax += 1

    d = []
    # line predictions
    for vg in range(0, vmax):
        for ve in range(0, vmax):
            for Jg in range(0, Jmax):
                for Je in range(0, Jmax):
                    
                    if (np.abs(Je - Jg) <= 1) and not (Je == 0):
                        eng = get_energy(Yg35, Ye35, vg, Jg, ve, Je)

                        hlp = [vg, Jg, ve, Je, eng, eng * 100.0 * c/1e12, eng * 100.0 * c/3e12]

                        d.append(hlp)

    f = open(save_filename, 'w')
    np.savetxt(f, d, delimiter = ',')
    f.close()

    print_matrix(d)

    # print order-reduced coefficients

    # (we  + wexe  * (v+1/2) + weye * (v+1/2)^2) * (v+1/2)
    # (we' + wexe' * (v+1/2)) * (v+1/2)
    # (we'') * (v+1/2)



    return 



