import numpy as np
from configparser import ConfigParser
import ast
import matplotlib.pyplot as plt

c = 299792458
kB = 1.3806482e-23
h_planck = 6.63607e-34
amu = 1.66e-27


def get_pop(J, B, T, N0 = 1):
    # returns rotational population distribution

    pop = N0 * (2*J+1) * np.exp(-B * h_planck * J*(J+1)/(kB*T))

    return pop / (np.sqrt(kB*T/(2*h_planck*B)) - 1.0/2.0)

def get_doppler(T, f0):
    # returns Doppler width for 27Al35Cl

    return np.sqrt(8*kB*T*np.log(2)/((35+27)*amu*c**2)) * f0


def simple_gauss(x, x0, A, w):
    return A * np.exp( -(x-x0)**2/(2*w**2) )


def av(arr, no_of_avg):
    # for 1D array
    if len(arr.shape)==1:
        hlp = np.zeros([int(arr.shape[0]/no_of_avg)])

        for k in range(len(hlp)):
            for m in range(no_of_avg):
                hlp[k] += arr[no_of_avg*k + m]

        return hlp/no_of_avg

    if len(arr.shape)==2:     

        # for 2D array
        hlp = np.zeros([int(arr.shape[0]/no_of_avg), arr.shape[1]])

        for k in range(len(hlp)):
            for m in range(no_of_avg):
                hlp[k] += arr[no_of_avg*k + m, :]

        return hlp/no_of_avg

def read_in_config(f):
    
    config = ConfigParser()
    config.read(f)

    sensor_ids = config.sections()
    # make dictionary out of config

    sensors = {}

    for s in sensor_ids:
        opts = config.options(s)
        
        sensors[s] = {}
        for o in opts:
            sensors[s][o] = config.get(s, o)

    return sensors

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n



def scale_coeff(M, mass1, mass2, scale = True):

    # scale or unscale Dunham coefficients with the reduced mass
    mu = (mass1 * mass2)/(mass1 + mass2)

    if scale == False:
        
        # return Y_kl = 1/mu^(-k/2 - l) * U_kl
        U = []
        for k in range(len(M)):
            hlp = []
            for l in range(len(M[k])):
                hlp.append(1.0/mu**(-k/2.0 - l) * M[k][l])
            U.append(hlp)
 
        return U
    else:

        # return Y_kl = mu^(-k/2 - l) * U_kl
        Y = []
        for k in range(len(M)):
            hlp = []
            for l in range(len(M[k])):
                hlp.append(mu**(-k/2.0 - l) * M[k][l])
            Y.append(hlp)
        
        return Y


def energy(Y, v, J):

    e = 0.0
    for k in range(len(Y)):
        for l in range(len(Y[k])):

            # E = Ykl(Be) * (v + 0.5)**k * ( J (J+1) )**l

            e += Y[k][l] * (v + 0.5)**k * ( J * (J + 1.0) )**l

    return e

def get_transitions(Yg, Ye, ve, vg, cnt_freq, df = 100e6, T = 1, Jmax = 1, no_of_points = 3000, real_amplitudes = True, pressure_broadening = 1.0):

    Jg_arr = np.arange(0, Jmax+1)    
    Je_arr = np.arange(0, Jmax+1)
    
    w = get_doppler(pressure_broadening * T, cnt_freq)
    
    nus = np.linspace(cnt_freq - df, cnt_freq + df, no_of_points)
    
    spectrum_P = np.zeros(len(nus))
    spectrum_Q = np.zeros(len(nus))
    spectrum_R = np.zeros(len(nus))
    
    f_P = []
    f_Q = []
    f_R = []
    
    for Jg in Jg_arr:
        for Je in Je_arr:
            
            # only dipole transitions
    
            eng = 100*c*(energy(Ye, ve, Je) - energy(Yg, vg, Jg))
    
            # apply population of ground states
            if real_amplitudes:
                A = get_pop(Jg, 100*c*Yg[0][1], T)
            else:
                A = 1.0
    
            if Je - Jg == -1: 
                spectrum_P += simple_gauss(nus, eng, A, w)
                f_P.append(eng)
                #plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'r' )
    
            if Je - Jg ==  0: 
                spectrum_Q += simple_gauss(nus, eng, A, w)
                f_Q.append(eng)
                #plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'g' )
    
            if Je - Jg == +1: 
                spectrum_R += simple_gauss(nus, eng, A, w)
                f_R.append(eng)
                #plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'b' )
    
    
    f_P = np.array(f_P)
    f_Q = np.array(f_Q)
    f_R = np.array(f_R)

    return (nus, [spectrum_P, spectrum_Q, spectrum_R], [f_P, f_Q, f_R])

def get_spec_lines(Yg, Ye, ve, vg, nus, line_type = 'Q', T = 1, Jmax = 1, real_amplitudes = True, cnt_freq = 100.0 * c * 38237.0):

    Jg_arr = np.arange(0, Jmax+1)
    Je_arr = np.arange(0, Jmax+1)
    
    w = get_doppler(T, cnt_freq)
    
    spectrum = np.zeros(len(nus))
    
    for Jg in Jg_arr:
        for Je in Je_arr:
            
            # only dipole transitions
    
            eng = 100*c*(energy(Ye, ve, Je) - energy(Yg, vg, Jg))
    
            # apply population of ground states
            if real_amplitudes:
                A = get_pop(Jg, 100*c*Yg[1][0], T)
            else:
                A = 1.0
    
            if Je - Jg == -1 and line_type == 'P': 
                spectrum += simple_gauss(nus, eng, A, w)
    
            if Je - Jg ==  0 and line_type == 'Q': 
                spectrum += simple_gauss(nus, eng, A, w)
    
            if Je - Jg == +1 and line_type == 'R': 
                spectrum += simple_gauss(nus, eng, A, w)
            
    

    return (nus, spectrum)



def plot_spectrum(nus, s, ve = 0, vg = 0, T = 0, cnt_freq = 0, style = '-', txt = '', abundance = 1.0):
    
    plt.plot( (nus - cnt_freq) / 1e9, abundance * s[0], 'r' + style, label = 'P' + txt)
    plt.plot( (nus - cnt_freq) / 1e9, abundance * s[1], 'g' + style, label = 'Q' + txt)
    plt.plot( (nus - cnt_freq) / 1e9, abundance * s[2], 'b' + style, label = 'R' + txt)

    plt.xlabel("Frequency (GHz) - {0:6.6f} THz".format(cnt_freq/1e12))

    plt.title("AlCl transitions @ T = {2:2} K for v = {0} -> v' = {1} (J->J')".format(vg, ve, T))

    plt.xlim(np.min(nus-cnt_freq)/1e9, np.max(nus-cnt_freq)/1e9)

    #plt.ylim(-0.1, 3.0)

    plt.legend()

def plot_transitions(f, cnt_freq = 0.0, cut = None, style = '-', txt = ''):
    
    if cut == None:
       cut = len(f[0])

    plt.plot((f[0][0:cut] - cnt_freq)/1e9, 'ro' + style, label = 'P' + txt)
    plt.plot((f[1][0:cut] - cnt_freq)/1e9, 'gx' + style, label = 'Q' + txt)
    plt.plot((f[2][0:cut] - cnt_freq)/1e9, 'bd' + style, label = 'R' + txt)

    plt.xlabel('Rotational number J')
    plt.ylabel("Frequency (GHz) - {0:6.6f} THz".format(cnt_freq/1e12))

    plt.legend(loc = 'upper right')



def combine_data(x_arr, y_arr, arr = [], sort = True):

    # concatenates data arrays

    if len(arr) == 0:
        arr = range(len(x_arr))

    x = []
    y = []
    for n in arr:

        x.extend(x_arr[n])
        y.extend(y_arr[n])

    x = np.array(x)
    y = np.array(y)

    if sort:

        ind = np.argsort(x)

        x = x[ind]
        y = y[ind]

    return (x, y)

def get_reduced_dunham():

    # ground state from Bernath
    Ug = [
            [0.0,3.71517408,-5.802142e-5],
            [1880.202,-9.575654e-2],
            [-32.012],
            [3.95186e-1],
            [-4.802e-3]
            ]
    
    # excited state
    Ue = [
           [38251.101190451634, 3.6400046656060443, 0.00], 
           [1784.3665982617565, 0.0*-0.2344361133562231],
           [-114.5238709407005], 
           [17.527497484431063], 
           [-6.316749292366754]
         ]

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


