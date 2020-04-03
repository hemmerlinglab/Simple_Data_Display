import numpy as np
import pickle
import matplotlib.pyplot as plt
from fit_func import *
from helper import *

# load data

arr = pickle.load( open( "data.pickle", "rb" ) )

f_arr = arr['f_arr']
sig_arr = arr['sig_arr']
data = arr['data']


c = 299792458

cnt_freq = 100.0 * c * 38237.00






# AlCl constants

(Ug, Ue) = get_reduced_dunham()

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)



Jmax = 30
T = 10 # Kelvin
#cnt_freq = 100.0*c*38237.0 # 1/cm
df = 100e9 # Hz
cl35_abund = 0.76
cl37_abund = 0.24

# vibrational states
ve = 0
vg = 0
(nus35, spectrum35, f_lines35) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, df, T, Jmax, pressure_broadening = 10.0)
(nus37, spectrum37, f_lines37) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, df, T, Jmax, pressure_broadening = 10.0)

plt.figure()


plot_spectrum(nus35, spectrum35, ve, vg, T, cnt_freq, style = '-', txt = ' Cl35', abundance = cl35_abund)
plot_spectrum(nus37, spectrum37, ve, vg, T, cnt_freq, style = '--', txt = ' Cl37', abundance = cl37_abund)


for n in range(8):

    cnt = np.mean(f_arr[n])

    plt.axvline((cnt - cnt_freq)/1e9, color = 'k')
    #plt.plot((f_arr[n] - cnt_freq)/1e9, sig_arr[n]*100.0, 'ko')



mu_35 = (massAl * massCl_35)/(massAl + massCl_35)
mu_37 = (massAl * massCl_37)/(massAl + massCl_37)

img = []



UBe_arr = np.linspace(3.7, 3.9, 25)


for k in range(len(Be_arr)):

    Ue[0][1] = UBe_arr[k]

    Ye35 = scale_coeff(Ue, massAl, massCl_35, scale = True)
    Ye37 = scale_coeff(Ue, massAl, massCl_37, scale = True)

    #(nus35, spectrum35, f_lines35) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, df, T, Jmax, pressure_broadening = 3.0)
    #(nus37, spectrum37, f_lines37) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, df, T, Jmax, pressure_broadening = 3.0)

    #img.append(cl35_abund * np.sum(spectrum35, axis = 0) + cl37_abund * np.sum(spectrum37, axis = 0))
    

    ve = 0
    vg = 0
    (nus35, spec35) = get_spec_lines(Yg35, Ye35, ve, vg, x, line_type = 'P', T = T, Jmax = Jmax, real_amplitudes = False)
    (nus37, spec37) = get_spec_lines(Yg37, Ye37, ve, vg, x, line_type = 'P', T = T, Jmax = Jmax, real_amplitudes = False)

    f = spec35 + spec37

    


plt.figure()
plt.pcolor((nus35 - cnt_freq)/1e9, Be_arr, img)
plt.clim(0,0.75)
plt.colorbar()

# add all data to the plot

for n in range(8):

    cnt = np.mean(f_arr[n])

    plt.axvline((cnt - cnt_freq)/1e9, color = 'k')
    #plt.plot((f_arr[n] - cnt_freq)/1e9, sig_arr[n]*100.0, 'ko')



plt.show()


