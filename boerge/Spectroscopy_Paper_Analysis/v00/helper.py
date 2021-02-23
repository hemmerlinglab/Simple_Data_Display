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


def energy(Y, v, J):

    e = 0.0
    for k in range(len(Y)):
        for l in range(len(Y[k])):

            e += Y[k][l] * (v + 0.5)**k * ( J * (J + 1.0) )**l

    return e

def get_transitions(Yg, Ye, ve, vg, cnt_freq, df = 100e9, T = 10, Jmax = 1, no_of_points = 3000, real_amplitudes = True, pressure_broadening = 1.0, isotope_abundance = 1.0):

    Jg_arr = np.arange(0, Jmax+1)    
    Je_arr = np.arange(1, Jmax+1)
    
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

    return (nus, [isotope_abundance * spectrum_P, isotope_abundance * spectrum_Q, isotope_abundance * spectrum_R], [f_P, f_Q, f_R])

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

    # these fit well for both transitions v=0 -> v=0 and v=1 -> v=1
    # they are far off from what Brian sees, though
    Ug = [[0, 3.697689639, -0.00004110220051, 0, 0, 0], [1880.20433, \
    0.04887253993, 0, 0, 0, 0], [-32.01271, 0, 0, 0, 0, 0], [0.0395499, \
    0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

    Ue = [[38255.63489, 3.718219212, 0.00273806453, 0, 0, 0], [1729.598167, \
    -0.06421374131, 0, 0, 0, 0], [60.74391688, 0, 0, 0, 0, 0], \
    [-180.1417929, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, \
    0]]

    ## less orders
    #Ug = [[0, 3.703932569, -0.00004535349545], [1880.20433, 0, 0], \
    #[-32.01271, 0, 0]]
    #Ue = [[38255.41096, 3.710086911, 0.002721267251], [1736.170635, 0, 0], [0, \
    #0, 0]]

    #Ug = [[0, 3.058994665, -0.00008158000278], [1880.20433, 5.063896456, 0], \
    #[-32.01271, 0, 0]]
    #Ue = [[38255.45438, 3.140849542, 0.002912864012], [1735.832448, \
    #4.464122387, 0], [0, 0, 0]]

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



def read_in_data(data, offset_avg = 10, datafolder = '/Users/boerge/Software/offline_data/', moving_avg_no = 0):
 
    basefolder = str(data['date'])

    basefilename = datafolder + basefolder + '/' + basefolder + '_'

    f_freqs = basefilename + str(data['time']) + '_set_points'
    f_act_freqs = basefilename + str(data['time']) + '_act_freqs'
    f_ch0   = basefilename + str(data['time']) + '_ch0_arr'
    f_ch1   = basefilename + str(data['time']) + '_ch1_arr'
    f_ch2   = basefilename + str(data['time']) + '_ch2_arr'
    f_ch3   = basefilename + str(data['time']) + '_ch3_arr'

    config_file = basefilename + str(data['time']) + '_conf'

    conf = read_in_config(config_file)

    print('Analyzing file ... ' + f_freqs)
    
    freqs = np.genfromtxt(f_freqs, delimiter=",")
    act_freqs = np.genfromtxt(f_act_freqs, delimiter=",")
    ch0 = np.genfromtxt(f_ch0, delimiter=",")
    ch1 = np.genfromtxt(f_ch1, delimiter=",")
    ch2 = np.genfromtxt(f_ch2, delimiter=",")
    ch3 = np.genfromtxt(f_ch3, delimiter=",")

    # get number of averages
    no_of_avg = int(len(freqs)/len(np.unique(freqs)))

    print('Found ' + str(no_of_avg) + ' averages.')

    # take the averages
    freqs = av(freqs, no_of_avg)
    act_avg_freqs = av(act_freqs, no_of_avg)
    ch0 = av(ch0, no_of_avg)
    ch1 = av(ch1, no_of_avg)
    ch2 = av(ch2, no_of_avg)
    ch3 = av(ch3, no_of_avg)

    # subtracting the DC offset
    offset_avg_points = offset_avg
    for k in range(ch0.shape[0]):
        ch0[k, :] = ch0[k, :] - np.mean(ch0[k, -offset_avg_points:-1])
        ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
        ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])
        ch3[k, :] = ch3[k, :] - np.mean(ch3[k, -offset_avg_points:-1])

    # get frequency scan interval in terms of absolute frequencies
    laser_offset = 3.0 * np.float(conf['offset_laser1']['val']) * 1e12

    act_freqs = 3.0 * act_freqs * 1e12
    act_avg_freqs = 3.0 * act_avg_freqs * 1e12
    freqs = 3.0 * freqs * 1e6
 
    # UV frequency = 3x IR frequency
    freqs = freqs + laser_offset

    t_steps = np.float(conf['step_size']['val'])*1e-6

    t_count = np.float(conf['scope_count']['val'])

    times = np.linspace(0, t_steps * t_count, np.int(t_count))

    times = times / 1e-3
   
    # apply moving time average 
    if moving_avg_no > 0:
        
        times = moving_average(times, n = moving_avg_no)

        ch0_avg = np.zeros([ch0.shape[0], len(times)])
        ch1_avg = np.zeros([ch1.shape[0], len(times)])
        ch2_avg = np.zeros([ch2.shape[0], len(times)])
        ch3_avg = np.zeros([ch3.shape[0], len(times)])
    
        for k in range(ch0.shape[0]):
            ch0_avg[k, :] = moving_average(ch0[k, :], n = moving_avg_no)
            ch1_avg[k, :] = moving_average(ch0[k, :], n = moving_avg_no)
            ch2_avg[k, :] = moving_average(ch0[k, :], n = moving_avg_no)
            ch3_avg[k, :] = moving_average(ch0[k, :], n = moving_avg_no)

        return (times, freqs, act_freqs, act_avg_freqs, [ch0_avg, ch1_avg, ch2_avg, ch3_avg], laser_offset)

    else:
        
        return (times, freqs, act_freqs, act_avg_freqs, [ch0, ch1, ch2, ch3], laser_offset)



def print_dunham():

    (Ug, Ue) = get_reduced_dunham()
    
    (Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)
    

    print('Dunham coefficients')
    print('-'*30)
    #for v in range(len(Yg35)):
    #    for l in range(len(Yg35[v])):

    print("Yg35_00 = {0:8.3f} THz".format(Yg35[0][0]*100*c/1e12))
    print("Yg35_01 = {0:8.3f} GHz".format(Yg35[0][1]*100*c/1e9))
    print("Yg35_02 = {0:8.3f} kHz".format(Yg35[0][2]*100*c/1e3))
    print("Yg35_10 = {0:8.3f} THz".format(Yg35[1][0]*100*c/1e12))
    print("Yg35_20 = {0:8.3f} GHz".format(Yg35[2][0]*100*c/1e9))
    print("Yg35_11 = {0:8.3f} MHz".format(Yg35[1][1]*100*c/1e6))
    print()
    print("Ye35_00 = {0:8.3f} THz".format(Ye35[0][0]*100*c/1e12))
    print("Ye35_01 = {0:8.3f} GHz".format(Ye35[0][1]*100*c/1e9))
    print("Ye35_02 = {0:8.3f} kHz".format(Ye35[0][2]*100*c/1e3))
    print("Ye35_10 = {0:8.3f} THz".format(Ye35[1][0]*100*c/1e12))
    print("Ye35_20 = {0:8.3f} GHz".format(Ye35[2][0]*100*c/1e9))
    print("Ye35_11 = {0:8.3f} MHz".format(Ye35[1][1]*100*c/1e6))




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
   


