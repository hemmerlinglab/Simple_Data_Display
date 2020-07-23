import numpy as np
from helper import *
import matplotlib.pyplot as plt
from fit_func import *



ind1 = 98
ind2 = 120 #105

data = [
        {'date':20200616, 'time':142554, 'ind1':ind1, 'ind2':ind2, 'color' : 'k', 'Jg' : 1, 'Je' : 1, 'iso' : 35, 'skip' : True}, # Q line
        {'date':20200616, 'time':144555, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 0, 'Je' : 1, 'iso' : 35},
        {'date':20200616, 'time':145926, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 0, 'Je' : 1, 'iso' : 37},
        {'date':20200616, 'time':150759, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 1, 'Je' : 2, 'iso' : 35},
        {'date':20200616, 'time':152334, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 1, 'Je' : 2, 'iso' : 37},
        {'date':20200616, 'time':153802, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 2, 'Je' : 3, 'iso' : 35},
        {'date':20200616, 'time':154644, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 2, 'Je' : 3, 'iso' : 37},
        {'date':20200616, 'time':155544, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 3, 'Je' : 4, 'iso' : 35},
        {'date':20200616, 'time':160401, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 3, 'Je' : 4, 'iso' : 37},
        {'date':20200616, 'time':161314, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 4, 'Je' : 5, 'iso' : 35},
        {'date':20200616, 'time':162146, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 4, 'Je' : 5, 'iso' : 37},
        {'date':20200616, 'time':162955, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 5, 'Je' : 6, 'iso' : 35},
        {'date':20200616, 'time':163757, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 5, 'Je' : 6, 'iso' : 37},
        {'date':20200616, 'time':170320, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 6, 'Je' : 7, 'iso' : 35},
        {'date':20200616, 'time':171158, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 6, 'Je' : 7, 'iso' : 37, 'skip' : True},
        {'date':20200616, 'time':172034, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 7, 'Je' : 8, 'iso' : 35, 'skip' : True},
        {'date':20200616, 'time':172939, 'ind1':ind1, 'ind2':ind2, 'color' : 'r', 'Jg' : 7, 'Je' : 8, 'iso' : 37, 'skip' : True},
        {'date':20200618, 'time':111956, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 2, 'Je' : 1, 'iso' : 37},
        {'date':20200618, 'time':120857, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 2, 'Je' : 1, 'iso' : 35},
        {'date':20200618, 'time':133509, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 3, 'Je' : 2, 'iso' : 37},
        {'date':20200618, 'time':141057, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 3, 'Je' : 2, 'iso' : 35},
        {'date':20200618, 'time':141959, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 4, 'Je' : 3, 'iso' : 37},
        {'date':20200618, 'time':142824, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 4, 'Je' : 3, 'iso' : 35},
        {'date':20200618, 'time':143538, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 5, 'Je' : 4, 'iso' : 37},
        {'date':20200618, 'time':144338, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 5, 'Je' : 4, 'iso' : 35},
        {'date':20200618, 'time':145042, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 6, 'Je' : 5, 'iso' : 37},
        {'date':20200618, 'time':145759, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 6, 'Je' : 5, 'iso' : 35},
        {'date':20200618, 'time':150527, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 7, 'Je' : 6, 'iso' : 37},
        {'date':20200618, 'time':151343, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 7, 'Je' : 6, 'iso' : 35},
        {'date':20200618, 'time':152045, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 8, 'Je' : 7, 'iso' : 37},
        {'date':20200618, 'time':153043, 'ind1':ind1, 'ind2':ind2, 'color' : 'b', 'Jg' : 8, 'Je' : 7, 'iso' : 35},
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
    n_mov = 6#4
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

#cnt_freq = 3*381.747875e12

for n in range(len(data)):
    plt.plot((f_arr[n] - cnt_freq)/1e9, sig_arr[n], 'o-', label = data[n]['time'], color = data[n]['color'])
    plt.text(np.mean((f_arr[n] - cnt_freq)/1e9), np.max(sig_arr[n]), str(n))

#plt.legend()

plt.xlabel("Frequency (GHz) + {0:3.6f} THz".format(cnt_freq/1e12))


############################# rough theory



(Ug, Ue) = get_reduced_dunham()

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)

ve = 0
vg = 0

(nus, [P35, Q35, R35], [f_P, f_Q, f_R]) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, df = 120e9, Jmax = 10, T = 4, real_amplitudes = True, isotope_abundance = 0.76)
(nus, [P37, Q37, R37], [f_P, f_Q, f_R]) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, df = 120e9, Jmax = 10, T = 4, real_amplitudes = True, isotope_abundance = 0.24)


nus = nus - cnt_freq
nus = nus/1e9

a = -0.006/4.0 * 0.5
o = -0.0001*0

plt.plot(nus, a * P35 + o, 'r')
plt.plot(nus, a * Q35 + o, 'k')
plt.plot(nus, a * R35 + o, 'b')

plt.plot(nus, a * P37 + o, 'r--')
plt.plot(nus, a * Q37 + o, 'k--')
plt.plot(nus, a * R37 + o, 'b--')


#a = 0.006/2
#o = 0
#
#plt.plot(nus, a * P35 + o, 'g-')
#plt.plot(nus, a * Q35 + o, 'g-')
#plt.plot(nus, a * R35 + o, 'g-')
#
#plt.plot(nus, a * P37 + o, 'r-')
#plt.plot(nus, a * Q37 + o, 'r-')
#plt.plot(nus, a * R37 + o, 'r-')











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



