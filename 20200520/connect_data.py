import numpy as np
from helper import *
import matplotlib.pyplot as plt
from fit_func import *



ind1 = 98
ind2 = 105

data = [
        {'date':20200520, 'time':132412, 'ind1':ind1, 'ind2':ind2, 'color' : 'k'},
        {'date':20200520, 'time':133143, 'ind1':ind1, 'ind2':ind2, 'color' : 'b'},
        {'date':20200520, 'time':141726, 'ind1':ind1, 'ind2':ind2, 'color' : 'b'},
        {'date':20200520, 'time':135733, 'ind1':ind1, 'ind2':ind2, 'color' : 'b'},
        {'date':20200520, 'time':140524, 'ind1':ind1, 'ind2':ind2, 'color' : 'b'},
        {'date':20200520, 'time':142617, 'ind1':ind1, 'ind2':ind2, 'color' : 'b'},
        {'date':20200520, 'time':143736, 'ind1':ind1, 'ind2':ind2, 'color' : 'b'},
        {'date':20200520, 'time':145052, 'ind1':ind1, 'ind2':ind2, 'color' : 'r'},
        {'date':20200520, 'time':145646, 'ind1':ind1, 'ind2':ind2, 'color' : 'k'},
        {'date':20200520, 'time':150858, 'ind1':ind1, 'ind2':ind2, 'color' : 'r'},
        {'date':20200520, 'time':152131, 'ind1':ind1, 'ind2':ind2, 'color' : 'r'},
        {'date':20200520, 'time':153414, 'ind1':ind1, 'ind2':ind2, 'color' : 'r'},
        {'date':20200520, 'time':154653, 'ind1':ind1, 'ind2':ind2, 'color' : 'r'},
                ]


f_arr = []
sig_arr = []
sig_fit_arr = []
hlp_arr = []

for n in range(len(data)):

    # read in data
    (times, freqs, ch_arr, laser_offset) = read_in_data(data[n], moving_avg_no = 0)

    #sig = np.mean(ch_arr[0], axis = 0)
    ## plot time traces
    #plt.figure()

    #plt.subplot(2,2,1)
    #plt.plot(times, sig)

    #plt.axvline(times[data[n]['ind1']])
    #plt.axvline(times[data[n]['ind2']])

    #plt.xlim(times[data[n]['ind1']] - 1.0, times[data[n]['ind2']] + 10)
    #plt.ylim(np.min(sig), 0)
    
    ind1 = data[n]['ind1']
    ind2 = data[n]['ind2']
    signal = -np.mean(ch_arr[0][:, ind1:ind2], axis = 1)

    # apply moving average
    n_mov = 4
    signal = moving_average(signal, n = n_mov)
    freqs = moving_average(freqs, n = n_mov)

    # append all signals to arrays
    f_arr.append(freqs)
    sig_arr.append(signal)
    #sig_fit_arr.append(hlp_fit)

f_arr = np.array(f_arr)
sig_arr = np.array(sig_arr)




plt.figure()

c = 299792458

#cnt_freq = 100.0 * c * 38237.00

cnt_freq = 1146.330906e12


for n in range(len(data)):
    plt.plot((f_arr[n] - cnt_freq)/1e9, sig_arr[n], 'o-', label = data[n]['time'], color = data[n]['color'])

#plt.legend()

plt.xlabel("Frequency (GHz) + {0:3.6f} THz".format(cnt_freq/1e12))



############################# rough theory



(Ug, Ue) = get_reduced_dunham()

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)

ve = 0
vg = 0
cnt_freq = 1146.330906e12

(nus, [P35, Q35, R35], [f_P, f_Q, f_R]) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, Jmax = 10, T = 4, real_amplitudes = False, isotope_abundance = 0.76)
(nus, [P37, Q37, R37], [f_P, f_Q, f_R]) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, Jmax = 10, T = 4, real_amplitudes = False, isotope_abundance = 0.24)

#print(energy(Yg35, 0, 0))
#print(energy(Ye35, 0, 0))

nus = nus - cnt_freq
nus = nus/1e9

a = -0.006/4.0 * 0.5
o = -0.001

plt.plot(nus, a * P35 + o, 'r')
plt.plot(nus, a * Q35 + o, 'k')
plt.plot(nus, a * R35 + o, 'b')

plt.plot(nus, a * P37 + o, 'r--')
plt.plot(nus, a * Q37 + o, 'k--')
plt.plot(nus, a * R37 + o, 'b--')


a = 0.006/2
o = 0

plt.plot(nus, a * P35 + o, 'g-')
plt.plot(nus, a * Q35 + o, 'g-')
plt.plot(nus, a * R35 + o, 'g-')

plt.plot(nus, a * P37 + o, 'g-')
plt.plot(nus, a * Q37 + o, 'g-')
plt.plot(nus, a * R37 + o, 'g-')











# save arrays

import pickle

with open('data.pickle', 'wb') as f:
    pickle.dump({'f_arr' : f_arr, 'sig_arr' : sig_arr, 'data' : data}, f)



# save as ascii file
x = []
y = []
for k in range(len(f_arr)):
    x.extend(f_arr[k])
    y.extend(sig_arr[k])


f = open('freqs.txt', 'w')
np.savetxt(f, x, delimiter = ',')
f.close()

f = open('signal.txt', 'w')
np.savetxt(f, y, delimiter = ',')
f.close()



print(len(x))
print(len(y))

plt.show()



