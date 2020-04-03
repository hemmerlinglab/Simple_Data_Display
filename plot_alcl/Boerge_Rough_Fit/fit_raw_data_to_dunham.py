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


# fit each data set with a gaussian to find center of lines (roughly)

line_pos = []

x_data = []
offset_free_data = []

######################################
# P lines (v -> v') = 0->0
######################################

#plt.figure()
for n in [0,1,2,3,4,5]:

    x = f_arr[n]
    y = sig_arr[n]

    x_mean = np.mean(x)
    x_min = np.min(x)
    x_max = np.max(x)
    y_min = np.min(y)
    y_max = np.max(y)

    params = Parameters()

    params.add('p0', value=y_min, min=-1.0, max=1.0, vary = True)
    params.add('p1', value=(y_max - y_min), min=0.0, max=3.0, vary = True)
    params.add('p2', value=x_mean, min=x_min-10e9, max=x_mean+10e9, vary = True)
    params.add('p3', value=(x_max-x_min)/10.0, min=0.001, max=(x_max-x_min), vary = True)
 
    (xf, yf, m) = fit_func(x, y, params, gauss)

    line_pos.append(m.params['p2'].value)

    x_data.append(x)
    offset_free_data.append((y - m.params['p0'])/m.params['p1'])

    #plt.subplot(3,3,n+1)
    #plt.plot((x - cnt_freq)/1e9, y, 'ro', label = data[n]['time'])
    #plt.plot((xf - cnt_freq)/1e9, yf, 'k-')

    #plt.xlabel("f (GHz) - {0:2.3f} THz".format(cnt_freq/1e12))
    #plt.legend()
    #plt.tight_layout()


#######################################################################
# Q line(s) (v -> v') = 0->0, need to be fitted to a sum of gaussians
#######################################################################


x = []
y = []
for n in [6,7]:

    x.extend(f_arr[n])
    y.extend(sig_arr[n])

x = np.array(x)
y = np.array(y)

x_mean = np.mean(x)
x_min = np.min(x)
x_max = np.max(x)
y_min = np.min(y)
y_max = np.max(y)

params = Parameters()

params.add('p0', value=y_min, min=-1.0, max=1.0, vary = True)
params.add('p3', value=(x_max-x_min)/10.0, min=0.001, max=(x_max-x_min), vary = True)

for k in range(4):
    # amplitudes
    params.add('p1' + str(k), value=(y_max - y_min), min=0.0, max=3.0, vary = True)

    # positions
    params.add('p2' + str(k), value=x_min + 1.0e9 * k, min=x_min-10e9, max=x_mean+10e9, vary = True)


(xf, yf, m) = fit_func(x, y, params, multi_gauss)

for k in range(4):
    line_pos.append(m.params['p2' + str(k)].value)

x_data.append(x)
offset_free_data.append( (y - m.params['p0']) / m.params['p10'] )

#plt.figure()
#plt.plot((x - cnt_freq)/1e9, y, 'ro')
#plt.plot((xf - cnt_freq)/1e9, yf, 'k-')

#plt.xlabel("f (GHz) - {0:2.3f} THz".format(cnt_freq/1e12))
#plt.tight_layout()




############################################################################
# Q line(s) (v -> v') = 0->0, need to be fitted to a sum of gaussians
#########################################################3#################

line_pos = np.sort(np.array(line_pos))

#plt.figure()
#plt.plot((line_pos - cnt_freq)/1e9, 'x')




line_pos = np.array(line_pos)

##################################################################################################
# get all fitted lines and subtract offset, then compare to simulated data by varying Be and Te
##################################################################################################

# take only unique values
uni_lines = line_pos[[0,1,2,5,6,8,9]]

gauss_data = []

ampl = 1.0
width = 0.5e9
gauss_x = cnt_freq + np.linspace(-100e9, 30e9, 500)

for k in range(len(uni_lines)):

    gauss_data.append(simple_gauss(gauss_x, uni_lines[k], ampl, width))

gauss_data = np.array(gauss_data)

x_data = np.array(x_data)
offset_free_data = np.array(offset_free_data)









#################################################################

# AlCl constants

(Ug, Ue) = get_reduced_dunham()

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)




###########################
# Fit the P-lines first
###########################


massAl = 26.98153841
massCl_35 = 34.96885269
massCl_37 = 36.96590258

mu35 = (massAl * massCl_35)/(massAl + massCl_35)
mu37 = (massAl * massCl_37)/(massAl + massCl_37)


# only P-lines
(x, y) = combine_data(x_data, offset_free_data, arr = [0,1,2,3,4,5])

T = 5
Jmax = 10

UBe_arr = 3.72 + np.linspace(-0.05, 0.05, 20)
UTe_arr = 38251.82 + np.linspace(-0.05, 0.05, 60)

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)

fit_min = []

for k in range(len(UBe_arr)):
    print(k)
    hlp = []
    for k2 in range(len(UTe_arr)):

       Ye35[0][1] = 1/mu35 * UBe_arr[k]
       Ye37[0][1] = 1/mu37 * UBe_arr[k]
       
       Ye35[0][0] = UTe_arr[k2]
       Ye37[0][0] = UTe_arr[k2]

       ve = 0
       vg = 0
       (nus35, spec35) = get_spec_lines(Yg35, Ye35, ve, vg, x, line_type = 'P', T = T, Jmax = Jmax, real_amplitudes = False)
       (nus37, spec37) = get_spec_lines(Yg37, Ye37, ve, vg, x, line_type = 'P', T = T, Jmax = Jmax, real_amplitudes = False)

       f = spec35 + spec37

       hlp.append(np.sum((f - y)**2))
    fit_min.append(hlp)


plt.figure()
plt.pcolor(UTe_arr - 38251.0, UBe_arr, fit_min)
plt.colorbar()


ind = np.where(fit_min == np.min(fit_min))

ind_Be = ind[0]
ind_Te = ind[1]


Ue[0][0] = UTe_arr[ind_Te]
Ue[0][1] = UBe_arr[ind_Be]

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)

xf = np.linspace(min(x), max(x), 1000)

ve = 0
vg = 0
(nus35, spec35) = get_spec_lines(Yg35, Ye35, ve, vg, xf, line_type = 'P', T = T, Jmax = Jmax, real_amplitudes = False)
(nus37, spec37) = get_spec_lines(Yg37, Ye37, ve, vg, xf, line_type = 'P', T = T, Jmax = Jmax, real_amplitudes = False)

f = spec35 + spec37

print(Ue)

plt.figure()
plt.plot((x - cnt_freq)/1e9, y, 'o')
plt.plot((xf - cnt_freq)/1e9, f, 'k')

plt.show()

asd
    

















params = Parameters()

params.add('p0', value = 38251.01, min = 38245.0, max = 38262.0, vary = True) # TeA
params.add('p1', value = 3.72, min = 0.2000, max = 5.000, vary = True) # BeA
params.add('p2', value = -0.23, min = -0.3, max = 1.1, vary = True) # Y11


#params.add('p0', value = 38477.8, min = 38470.0, max = 385002.0, vary = True) # TeA
#params.add('p1', value = 3.74556, min = 0.2000, max = 7.000, vary = True) # BeA

#params.add('p2', value = -0.003944, min = -0.0045, max = -0.0030, vary = True) # -De

(xf, yf, m) = fit_func(x, y, params, P_lines)

(yinit) = P_lines(xf, params)

print_results(m)


plt.figure()
plt.plot((x - cnt_freq)/1e9, y, 'ro')
plt.plot((xf - cnt_freq)/1e9, yf, 'k')
plt.plot((xf - cnt_freq)/1e9, yinit, 'b--')

plt.show()

asd




Jmax = 20
T = 20 # Kelvin
#cnt_freq = 100.0*c*38237.0 # 1/cm
df = 100e9 # Hz
cl35_abund = 0.76
cl37_abund = 0.24

# vibrational states
ve = 0
vg = 0
(nus35, spectrum35, f_lines35) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, df, T, Jmax, real_amplitudes = False)
(nus37, spectrum37, f_lines37) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, df, T, Jmax, real_amplitudes = False)

plt.figure()
plot_spectrum(nus35, spectrum35, ve, vg, T, cnt_freq, style = '-', txt = ' Cl35', abundance = cl35_abund)
plot_spectrum(nus37, spectrum37, ve, vg, T, cnt_freq, style = '--', txt = ' Cl37', abundance = cl37_abund)


# add all data to the plot
#for n in range(8):
#    plt.plot((f_arr[n] - cnt_freq)/1e9, sig_arr[n]*100.0, 'ko')

for n in range(len(x_data)):
    plt.plot((x_data[n] - cnt_freq)/1e9, offset_free_data[n], 'ko')

#plt.plot((gauss_x - cnt_freq)/1e9, np.sum(gauss_data, axis = 0))



plt.show()


asd


#plt.plot((nus35 - cnt_freq)/1e9, cl35_abund * np.sum(spectrum35, axis = 0) + cl37_abund*np.sum(spectrum37, axis = 0), 'k-')

plt.figure()
plot_transitions(f_lines35, cnt_freq = cnt_freq, cut = 12, style = '-', txt = ' Cl35')
plot_transitions(f_lines37, cnt_freq = cnt_freq, cut = 12, style = '--', txt = ' Cl37')


#plt.figure()

# calculate with varying value for the rotational constant


(nus35, spectrum35, f_lines35) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, df, T, Jmax)
(nus37, spectrum37, f_lines37) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, df, T, Jmax)


plt.show()




