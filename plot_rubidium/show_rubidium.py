import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys


c = 299792458


# hlp is helper variable



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



my_today = datetime.datetime.today()

#datafolder = '/Users/boerge/software/data/molecule_computer/'
#datafolder = '/home/molecules/software/data/'

datafolder = 'data/'

basefolder = '20191126'


time_stamp = '152843' #sys.argv[1] #'145643'

#basefilename = datafolder + basefolder + '/' + basefolder + '_'
basefilename = datafolder + '/' + basefolder + '_'


f_freqs = basefilename + time_stamp + '_set_points'
freqs = np.genfromtxt(f_freqs, delimiter=",")
# get number of averages
no_of_avg = int(len(freqs)/len(np.unique(freqs)))
print('Found ' + str(no_of_avg) + ' averages.')

for k in range(3):
    exec('f_ch' + str(k) + ' = "' + basefilename + time_stamp + '_ch' + str(k) + '_arr"')
    exec('ch' + str(k) + ' = np.genfromtxt(f_ch' + str(k) + ', delimiter=",")')
    exec('ch' + str(k) + ' = av(ch' + str(k) + ', no_of_avg)')

freqs = av(freqs, no_of_avg)

ch0_mean = np.mean(ch0, axis = 1)
ch1_mean = np.mean(ch1, axis = 1)
ch2_mean = np.mean(ch2, axis = 1)

ch0_mean -= ch0_mean[0]
ch1_mean -= ch1_mean[0]
ch2_mean -= ch2_mean[0]

#ch0_mean -= ch2_mean
#ch1_mean -= ch2_mean

diff_sig = ch0_mean - ch1_mean


delay_in_for_loop = 100e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3

   
Rb_85_D2_abs_freqs = 384.230406373
Rb_87_D2_abs_freqs = 384.230484468

v_abs = []
mytext = []
ind = []

# including hyperfine structure
v_abs.append(Rb_85_D2_abs_freqs + (-1.264888e9 + 100.205e6)/1e12)
mytext.append('Rb-85, 3->4')
v_abs.append(Rb_85_D2_abs_freqs + (-1.264888e9 - 20.435e6)/1e12)
mytext.append('Rb-85, 3->3')
v_abs.append(Rb_85_D2_abs_freqs + (-1.264888e9 - 83.835e6)/1e12)
mytext.append('Rb-85, 3->2')

v_abs.append(Rb_85_D2_abs_freqs + (+1.770843e9 - 20.435e6)/1e12)
mytext.append('Rb-85, 2->3')
v_abs.append(Rb_85_D2_abs_freqs + (+1.770843e9 - 83.835e6)/1e12)
mytext.append('Rb-85, 2->2')
v_abs.append(Rb_85_D2_abs_freqs + (+1.770843e9 - 113.208e6)/1e12)
mytext.append('Rb-85, 2->1')

v_abs.append(Rb_87_D2_abs_freqs + (-2.563005e9 + 193.7407e6)/1e12)
mytext.append('Rb-87, 2->3')
v_abs.append(Rb_87_D2_abs_freqs + (-2.563005e9 - 72.9112e6)/1e12)
mytext.append('Rb-87, 2->2')
v_abs.append(Rb_87_D2_abs_freqs + (-2.563005e9 - 229.8518e6)/1e12)
mytext.append('Rb-87, 2->1')

v_abs.append(Rb_87_D2_abs_freqs + (+4.271676e9 - 72.911e6)/1e12)
mytext.append('Rb-87, 1->2')
v_abs.append(Rb_87_D2_abs_freqs + (+4.271676e9 - 229.8518e6)/1e12)
mytext.append('Rb-87, 1->1')


setpoint_offset = 384.230
set_abs_freqs = freqs*1e6/1e12 + setpoint_offset

for k in range(len(v_abs)):
    ind.append( (v_abs[k] - setpoint_offset)*1.0e6 )

mymax = np.max(diff_sig)

plt.figure()
plt.subplot(2,1,1)
plt.plot(freqs, ch0_mean)
plt.plot(freqs, ch1_mean)
plt.plot(freqs, ch2_mean)

dv = 1.0/len(v_abs)
for k in range(len(v_abs)):
    plt.axvline(ind[k], ls = '--', color = 'r')
    plt.text(ind[k], -dv * k, mytext[k])

plt.subplot(2,1,2)
plt.plot(freqs, diff_sig)
for k in range(len(v_abs)):
    plt.axvline(ind[k], ls = '--', color = 'r')
    plt.text(ind[k], -dv * k, mytext[k])

plt.xlabel('Freqs (MHz) + ' + str(setpoint_offset) + ' THz')

plt.show()




