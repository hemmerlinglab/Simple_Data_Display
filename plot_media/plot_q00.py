import numpy as np
import matplotlib.pyplot as plt
from fit_dunham import energy

c = 299792458
kB = 1.3806482e-23
h_planck = 6.63607e-34
amu = 1.66e-27

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

def get_pop(J, B, T, N0 = 1):

    pop = N0 * (2*J+1) * np.exp(-B * h_planck * c * J*(J+1)/(kB*T))

    return pop / (np.sqrt(kB*T/(2*h_planck*c*B)) - 1/2)

def get_doppler(T, f0):

    return np.sqrt(8*kB*T*np.log(2)/((35+27)*amu*c**2)) * f0


def gauss(x, x0, A, w):
    return A * np.exp( -(x-x0)**2/(2*w**2) )


def get_transitions(Yg, Ye, ve, vg, cnt_freqs, df = 100e6, T = 1, Jmax = 1):

    Jg_arr = np.arange(0, Jmax+1)    
    Je_arr = np.arange(0, Jmax+1)
    
    w = get_doppler(T, 100*c*cnt_freq)
    
    nus = np.linspace(100*c*(cnt_freq) - df, 100*c*(cnt_freq) + df, 2000)
    
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
            A = get_pop(Jg, Yg[0][1], T)
    
            if Je - Jg == -1: 
                spectrum_P += gauss(nus, eng, A, w)
                f_P.append(eng)
                #plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'r' )
    
            if Je - Jg ==  0: 
                spectrum_Q += gauss(nus, eng, A, w)
                f_Q.append(eng)
                #plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'g' )
    
            if Je - Jg == +1: 
                spectrum_R += gauss(nus, eng, A, w)
                f_R.append(eng)
                #plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'b' )
    
    
    f_P = np.array(f_P) - 100*c*cnt_freq
    f_Q = np.array(f_Q) - 100*c*cnt_freq
    f_R = np.array(f_R) - 100*c*cnt_freq

    nus = nus - 100*c*cnt_freqs

    return (nus, [spectrum_P, spectrum_Q, spectrum_R], [f_P, f_Q, f_R])


def plot_spectrum(nus, s, ve = 0, vg = 0, T = 0, cnt_freq = 0, style = '-', txt = '', abundance = 1.0):
    
    #plt.plot( nus / 1e9, abundance * s[0], 'r' + style, label = 'P' + txt)
    plt.plot( nus / 1e9, abundance * s[1], 'g' + style, label = 'Q' + txt)
    #plt.plot( nus / 1e9, abundance * s[2], 'b' + style, label = 'R' + txt)

    plt.xlabel("Frequency (GHz) - {0:6.6f} THz".format(100*c*cnt_freq/1e12))

    plt.title("AlCl transitions @ T = {2:2} K for v = {0} -> v' = {1} (J->J')".format(vg, ve, T))

    plt.xlim(np.min(nus)/1e9, np.max(nus)/1e9)

    plt.ylim(-0.1, 3.0)

    plt.legend()

def plot_transitions(f, cut = None, style = '-', txt = ''):
    
    if cut == None:
       cut = len(f[0])

    plt.plot(f[0][0:cut]/1e9, 'ro' + style, label = 'P' + txt)
    plt.plot(f[1][0:cut]/1e9, 'gx' + style, label = 'Q' + txt)
    plt.plot(f[2][0:cut]/1e9, 'bd' + style, label = 'R' + txt)

    plt.xlabel('Rotational number J')
    plt.ylabel("Frequency (GHz) - {0:6.6f} THz".format(100*c*cnt_freq/1e12))

    plt.legend(loc = 'upper right')



def scale_coeff(Y, mass1, mass2, scale = True):

    # scale coefficients with the reduced mass
    mu = (mass1 * mass2)/(mass1 + mass2)

    nvib = len(Y[0])

    if scale == False:

        # return U_kl = 1/mu^(-k/2 - l) * Y_kl
        U = []
        hlp = []
        for k in range(nvib):
            l = 0
            hlp.append( 1/(mu**(-k/2 - l)) * Y[0][k] )
        U.append(hlp)

        hlp = []
        for k in [0,1]:
            l = 1
            hlp.append( 1/(mu**(-k/2 - l)) * Y[1][k] )
        U.append(hlp)

        return U
    else:

        # return Y_kl = mu^(-k/2 - l) * U_kl

        U = Y.copy()
        Y = []
        hlp = []
        for k in range(nvib):
            l = 0
            hlp.append( (mu**(-k/2 - l)) * U[0][k] )
        Y.append(hlp)

        hlp = []
        for k in [0,1]:
            l = 1
            hlp.append( (mu**(-k/2 - l)) * U[1][k] )
        Y.append(hlp)

        return Y




###########################################################################################################################################################

# AlCl constants

Yg = [[0.0, 482.0794801498224, -2.0285504669476797, 4.2244382524958823e-07, -1.4200906251997294e-05], [0.23086415576356475, -0.0015452837298277023]]
Ye = [[38250.600844308836, 457.2272074915964, -7.519552217461303, 0.29489258632225485, -0.02723237768418022], [0.23713820022739715, -0.003944286505039427]]

# from Bernath
Yg = [[0.0, 481.774655, -2.101, 6.638e-3, -2.025e-5],[0.2439300, -1.611e-3]]
Ye = [[38250.600844308836, 457.2272074915964, -7.519552217461303, 0.29489258632225485, -0.02723237768418022], [0.23713820022739715, -0.003944286505039427]]

massAl = 26.98153841
massCl_35 = 34.96885269
massCl_37 = 36.96590258

Ug = scale_coeff(Yg, massAl, massCl_35, scale = False)
Ue = scale_coeff(Ye, massAl, massCl_35, scale = False)

Yg35 = scale_coeff(Ug, massAl, massCl_35, scale = True)
Ye35 = scale_coeff(Ue, massAl, massCl_35, scale = True)

Yg37 = scale_coeff(Ug, massAl, massCl_37, scale = True)
Ye37 = scale_coeff(Ue, massAl, massCl_37, scale = True)






Jmax = 11
T = 40 # Kelvin
cnt_freq = 38237.0 # 1/cm
df = 10e9 # HZ

# vibrational states
ve = 0
vg = 0
(nus35, spectrum35, f_lines35) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, df, T, Jmax)
(nus37, spectrum37, f_lines37) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, df, T, Jmax)

#plt.figure()
#plot_spectrum(nus35, spectrum35, ve, vg, T, cnt_freq, style = '-', txt = ' Cl35', abundance = 0.76)
#plot_spectrum(nus37, spectrum37, ve, vg, T, cnt_freq, style = '--', txt = ' Cl37', abundance = 0.24)



# fit q00 line



#datafolder = '/Users/boerge/software/data/molecule_computer/'
#datafolder = '/Users/boerge/Software/data/Molecules/Data/molecule_computer/'
datafolder = '/Users/boerge/tmp/'

basefolder = '20200227'

time_stamp = '111627'

#basefilename = datafolder + basefolder + '/' + basefolder + '_'
basefilename = datafolder + basefolder + '_'

f_freqs = basefilename + time_stamp + '_set_points'
f_ch1 = basefilename + time_stamp + '_ch0_arr'
f_ch2 = basefilename + time_stamp + '_ch1_arr'
f_ch4 = basefilename + time_stamp + '_ch3_arr'

freqs = np.genfromtxt(f_freqs, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")
ch2 = np.genfromtxt(f_ch2, delimiter=",")
ch4 = np.genfromtxt(f_ch4, delimiter=",")


# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
    ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
    ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])
    ch4[k, :] = ch4[k, :] - np.mean(ch4[k, -offset_avg_points:-1])




# get number of averages
no_of_avg = int(len(freqs)/len(np.unique(freqs)))

print('Found ' + str(no_of_avg) + ' averages.')

freqs = av(freqs, no_of_avg)
ch1 = av(ch1, no_of_avg)
ch2 = av(ch2, no_of_avg)
ch4 = av(ch4, no_of_avg)




time_stamp = '121347'

basefilename = datafolder + basefolder + '_'

f_freqs = basefilename + time_stamp + '_set_points'
f_ch1 = basefilename + time_stamp + '_ch0_arr'
f_ch2 = basefilename + time_stamp + '_ch1_arr'
f_ch4 = basefilename + time_stamp + '_ch3_arr'

freqs2 = np.genfromtxt(f_freqs, delimiter=",")
ch12 = np.genfromtxt(f_ch1, delimiter=",")
ch22 = np.genfromtxt(f_ch2, delimiter=",")
ch42 = np.genfromtxt(f_ch4, delimiter=",")


# subtracting the DC offset
offset_avg_points = 5
for k in range(ch12.shape[0]):
    ch12[k, :] = ch12[k, :] - np.mean(ch12[k, -offset_avg_points:-1])
    ch22[k, :] = ch22[k, :] - np.mean(ch22[k, -offset_avg_points:-1])
    ch42[k, :] = ch42[k, :] - np.mean(ch42[k, -offset_avg_points:-1])



# get number of averages
no_of_avg = int(len(freqs2)/len(np.unique(freqs2)))

print('Found ' + str(no_of_avg) + ' averages.')

freqs2 = av(freqs2, no_of_avg)
ch12 = av(ch12, no_of_avg)
ch22 = av(ch22, no_of_avg)
ch42 = av(ch42, no_of_avg)






freqs = np.append(freqs, freqs2)
sig = np.vstack((ch1, ch12))





#plt.figure()
#
#plt.plot(ch12[:, 0])
#plt.plot(ch22[:, 0])
#plt.plot(ch42[:, 0])
#
#plt.figure()
#
#plt.plot(ch12[0, :])
#plt.plot(ch22[0, :])
#plt.plot(ch42[0, :])



#nus = freqs*1.5


#delay_in_for_loop = 100e-6
delay_in_for_loop = 50e-6

no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3



#a1 = 3
#d1 = 20
#
#for d1 in range(1, 50, 10):
#    plt.figure()
#
#    for a1 in range(0,10,1):
#        plt.plot(freqs, np.mean(sig[:, 177+a1:177+a1+d1], axis = 1))
#
#plt.figure()


a1 = 1 #1
d1 = 10 #50
mean_sig = np.mean(sig[:, 177+a1:177+a1+d1], axis = 1)
mean_sig = mean_sig/np.min(mean_sig)




# fit the line
x_data = 3 * freqs
y_data = mean_sig

from fit_q00 import *
from fit_dl import *

(x_fit, y_fit, result) = fit_q00(x_data, y_data, Yg35, Ye35)

(x_fit2, y_fit2, result2) = fit_dl(x_data, y_data)


print(result.params)
print(result2.params)



plt.figure()
plt.plot(x_data/1e3, y_data, 'o')

plt.plot(x_fit2/1e3, y_fit2, 'k-')



plt.xlim(min(x_data/1e3), max(x_data/1e3))


plt.xlabel('Frequency (GHz)')
plt.ylabel('Absorption signal (a.u.)')


plt.show()




