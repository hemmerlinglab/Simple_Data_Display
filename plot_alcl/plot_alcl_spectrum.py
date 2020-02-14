import numpy as np
import matplotlib.pyplot as plt
from fit_dunham import energy

c = 299792458
kB = 1.3806482e-23
h_planck = 6.63607e-34
amu = 1.66e-27

def get_pop(J, B, T, N0 = 1):

    pop = N0 * (2*J+1) * np.exp(-B * h_planck * c * J*(J+1)/(kB*T))

    return pop / (np.sqrt(kB*T/(2*h_planck*c*B)) - 1/2)

def get_doppler(T, f0):

    return np.sqrt(8*kB*T*np.log(2)/((35+27)*amu*c**2)) * f0


def gauss(x, x0, A, w):
    return A * np.exp( -(x-x0)**2/(2*w**2) )


# AlCl constants

Yg = [[0.0, 482.0794801498224, -2.0285504669476797, 4.2244382524958823e-07, -1.4200906251997294e-05], [0.23086415576356475, -0.0015452837298277023]]
Ye = [[38250.600844308836, 457.2272074915964, -7.519552217461303, 0.29489258632225485, -0.02723237768418022], [0.23713820022739715, -0.003944286505039427]]




# calculate various lines



plt.figure()

Jmax = 30

T = 50 # Kelvin

cnt_freq = 38237.0

df = 100e9


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

ve = 0
vg = 0
for Jg in Jg_arr:
    for Je in Je_arr:
        
        # only dipole transitions

        eng = 100*c*(energy(Ye, ve, Je) - energy(Yg, vg, Jg))

        #if np.abs(Je - Jg) < 2: plt.text( (eng - 100*c*cnt_freq) / 1e9, 1.5, "{0}->{1}".format(Jg, Je))
        
        #if np.abs(Je - Jg) < 2: plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'k' )

        # apply population of ground states
        A = get_pop(Jg, Yg[0][1], T)

        if Je - Jg == -1: 
            spectrum_P += gauss(nus, eng, A, w)
            f_P.append(eng)
            plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'r' )

        if Je - Jg ==  0: 
            spectrum_Q += gauss(nus, eng, A, w)
            f_Q.append(eng)
            plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'g' )

        if Je - Jg == +1: 
            spectrum_R += gauss(nus, eng, A, w)
            f_R.append(eng)
            plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'b' )


f_P = np.array(f_P) - 100*c*cnt_freq
f_Q = np.array(f_Q) - 100*c*cnt_freq
f_R = np.array(f_R) - 100*c*cnt_freq


plt.plot( (nus - 100*c*cnt_freq) / 1e9, spectrum_P, 'r', label = 'P trans.')
plt.plot( (nus - 100*c*cnt_freq) / 1e9, spectrum_Q, 'g', label = 'Q trans.')
plt.plot( (nus - 100*c*cnt_freq) / 1e9, spectrum_R, 'b', label = 'R trans.')

plt.xlabel("Frequency (GHz) - {0:6.6f} THz".format(100*c*cnt_freq/1e12))

plt.title("AlCl transitions @ T = {2:2} K for v = {0} -> v' = {1} (J->J')".format(vg, ve, T))

plt.xlim((np.min(nus) - 100*c*cnt_freq)/1e9, (np.max(nus) - 100*c*cnt_freq)/1e9)

plt.ylim(-0.1, 3.0)

plt.legend()



# add second isotope
# using this paper: https://reader.elsevier.com/reader/sd/pii/0022285280901526?token=3D9E469371ED6A77EE19559548A05CF6AA45FEE434EF04A3DB0E380F9ADD8847CD9DD405AA148F56E3CD0AF6613F3404
# Journal of Molecular Spectroscopy 80, 411-421 (1980)


plt.figure()
plt.plot(f_P/1e9, 'ro-')
plt.plot(f_Q/1e9, 'gx-')
plt.plot(f_R/1e9, 'bd-')

plt.xlabel('Rotational number J')
plt.ylabel("Frequency (GHz) - {0:6.6f} THz".format(100*c*cnt_freq/1e12))




plt.figure()


cut = 5
plt.plot(f_P[0:cut]/1e9, 'ro-')
plt.plot(f_Q[0:cut]/1e9, 'gx-')
plt.plot(f_R[0:cut]/1e9, 'bd-')


exp_freqs = np.array([382.08140	 	 ,
382.08615	,
382.09095	,
382.09340	,
382.09575	,
382.09810	
])


exp_freqs = 3*exp_freqs*1e12

exp_freqs = (exp_freqs - 100*c*cnt_freq)

exp_freqs = exp_freqs/1e9

plt.plot(4, exp_freqs[0], 'kx', markersize = 10)
plt.plot(3, exp_freqs[1], 'kx', markersize = 10)
plt.plot(2, exp_freqs[2], 'kx', markersize = 10)
plt.plot(1, exp_freqs[3], 'kx', markersize = 10)
plt.plot(1, exp_freqs[4], 'kx', markersize = 10)
plt.plot(0, exp_freqs[5], 'kx', markersize = 10)


plt.xlabel('Rotational number J')
plt.ylabel("Frequency (GHz) - {0:6.6f} THz".format(100*c*cnt_freq/1e12))












plt.show()




