import numpy as np
import pickle
import matplotlib.pyplot as plt
from fit_func import *
from helper import *


c = 299792458

cnt_freq = 100 * c * 38237.0

arr = pickle.load( open( "line_centers.pickle", "rb" ) )

    
line_pos = arr['line_pos']
unique_lines = arr['unique_lines']
x_data = arr['x_data']
offset_free_data = arr['offset_free_data']


# take only the P-lines

#x_data = x_data[[0,1,2,3,4,5]]
#unique_lines = unique_lines[[0,1,2,3,4,5]]

#plt.figure()
#
#for n in range(len(x_data)):
#    plt.plot((x_data[n]-cnt_freq)/1e9, offset_free_data[n])
#
#print((unique_lines-cnt_freq)/1e9)
#plt.show()



massAl = 26.98153841
massCl_35 = 34.96885269
massCl_37 = 36.96590258

mu35 = (massAl * massCl_35)/(massAl + massCl_35)
mu37 = (massAl * massCl_37)/(massAl + massCl_37)





#################################################################

# AlCl constants

(Ug, Ue) = get_reduced_dunham()

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)


ve = 0
vg = 0
Jmax = 5

Jarr = np.arange(0, Jmax+1)


exp_lines = np.sort(unique_lines)



UBe_arr = 3.72 + np.linspace(-0.05, 0.05, 20)
UBe_arr = 3.72 + np.linspace(-0.15, 0.05, 40)
UBe_arr = 3.72 + np.linspace(-0.25, 0.05, 40)
UBe_arr = np.linspace(3.4, 3.6, 25)
UBe_arr = np.linspace(2.5, 4.0, 125)


cnt_test = 38251.82

UTe_arr = cnt_test + np.linspace(-0.05, 0.05, 30)
UTe_arr = cnt_test + np.linspace(-0.1, 0.1, 60)
UTe_arr = cnt_test + np.linspace(-0.2, 0.2, 60)
UTe_arr = cnt_test + np.linspace(-0.5, -0.2, 60)
UTe_arr = cnt_test + np.linspace(-0.35, -0.2, 25)

UTe_arr = cnt_test + np.linspace(-1.5, +1.5, 75)


fit_min = []

for k in range(len(UBe_arr)):
    print(k)
    hlp = []
    for k2 in range(len(UTe_arr)):

        eng = []
        # first the P-lines
        for Je in Jarr:
            
            # P transition
            Jg = Je + 1
        
            Ye35[0][1] = 1/mu35 * UBe_arr[k]
            Ye37[0][1] = 1/mu37 * UBe_arr[k]
               
            Ye35[0][0] = UTe_arr[k2]
            Ye37[0][0] = UTe_arr[k2]
        
            eng.append(100*c*(energy(Ye35, ve, Je) - energy(Yg35, vg, Jg)))
            eng.append(100*c*(energy(Ye37, ve, Je) - energy(Yg37, vg, Jg)))

        # then the Q line
        eng.append(100*c*(energy(Ye35, ve, 0) - energy(Yg35, vg, 0)))
        
        eng = np.sort(np.array(eng))

        # find for each line the closes in the array and build the subtraction
        sum_hlp = 0.0
        for n in range(len(exp_lines)):
            min_arr = np.abs(eng - exp_lines[n])

            ind = np.where( min_arr == np.min(min_arr) )

            sum_hlp += np.abs(eng[ind] - exp_lines[n])[0]/1e9

        hlp.append(sum_hlp)

    fit_min.append(hlp)

plt.figure()
plt.pcolor(UTe_arr - cnt_test, UBe_arr, fit_min)
plt.colorbar()









# plot energies with fitted parameters
# find minimum

ind = np.where(fit_min == np.min(np.min(fit_min)))

min_Be = UBe_arr[ind[0]]
min_Te = UTe_arr[ind[1]]

print(min_Be)
print(min_Te)

Ye35[0][1] = 1/mu35 * min_Be
Ye37[0][1] = 1/mu37 * min_Be
               
Ye35[0][0] = min_Te
Ye37[0][0] = min_Te


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
for n in range(len(x_data)):
    plt.plot((x_data[n] - cnt_freq)/1e9, offset_free_data[n], 'ko')




plt.show()



