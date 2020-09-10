import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
from fit_line import *

c = 299792458


def subtract_noise(data_array, noise_array):
    scls = []
    r = 80
    for i in range(r):
        scls.append(data_array[1,i]/noise_array[1,i])
    scl = np.mean(scls)
    noiseless = data_array - scl*noise_array
    plt.figure()
    plt.plot(np.mean(data_array[:, 1:r], axis = 1))
    plt.plot(np.mean(scl*noise_array[:, 1:r], axis = 1))

    return noiseless


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

#datafolder = '/Users/johnr/Desktop/'

datafolder = '/home/molecules/software/data/'



basefolder = '20200903'

basefilename = datafolder + basefolder + '/' + basefolder + '_'

if len(sys.argv)>1:
    time_stamp = sys.argv[1]
else:
    # get latest time stamp
    all_files = np.sort(glob.glob(basefilename + "*"))
    #print(all_files)
    time_stamp = all_files[-1].split('_')[1]

f_freqs = basefilename + time_stamp + '_set_points'
f_ch1 = basefilename + time_stamp + '_ch0_arr'
f_ch2 = basefilename + time_stamp + '_ch1_arr'
f_ch4 = basefilename + time_stamp + '_ch3_arr'

freqs = np.genfromtxt(f_freqs, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")
ch2 = np.genfromtxt(f_ch2, delimiter=",")
ch4 = np.genfromtxt(f_ch4, delimiter=",")

# get number of averages
no_of_avg = int(len(freqs)/len(np.unique(freqs)))

print('Found ' + str(no_of_avg) + ' averages.')

freqs = av(freqs, no_of_avg)
ch1 = av(ch1, no_of_avg)
ch2 = av(ch2, no_of_avg)
ch4 = av(ch4, no_of_avg)



nus = freqs*3.0

print(nus)
#delay_in_for_loop = 100e-6
delay_in_for_loop = 100e-6

no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3



cut_time1 = 10.00
cut_time2 = 11.00

ch1_start = np.where( np.abs(times - cut_time1) < 0.1 )[0][0]
ch1_end = np.where( np.abs(times - cut_time2) < 0.1 )[0][0]

# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
    ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
    ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])
    ch4[k, :] = ch4[k, :] - np.mean(ch4[k, -offset_avg_points:-1])




signal = np.mean(ch1[:, ch1_start:ch1_end], axis = 1)
#signal = 1-signal/np.min(signal)


mycut = 10
mycut2 = 300

abs_signal = np.mean(ch1[15:25, mycut:mycut2], axis = 0)
abs_signal = np.mean(ch1[:, mycut:mycut2], axis = 0)
abs_signal = 1.0 - abs_signal/np.min(abs_signal)


abs_times = times[mycut:mycut2]


myfontsize = 15

plt.figure()

(x, y, result) = fit_line(nus, signal) 

plt.subplot(2,1,1)
plt.plot(nus, signal, 'o-', markersize = 7)
plt.plot(x, y, 'r-', lw = 3)

plt.xlabel('Frequency (MHz)', fontsize = myfontsize)
plt.ylabel('Absorption signal (a.u.)', fontsize = myfontsize)

plt.xticks(fontsize = myfontsize)
plt.yticks(fontsize = myfontsize)

plt.tight_layout()


plt.subplot(2,1,2)
plt.plot(abs_times, abs_signal * 100, '-', lw = 3)

plt.xlabel('Times (ms)', fontsize = myfontsize)
plt.ylabel('Absorption signal (a.u.)', fontsize = myfontsize)


plt.ylim(-0.05 * 100,1.15 * 100)
plt.xlim(min(abs_times), max(abs_times))

plt.xticks(fontsize = myfontsize)
plt.yticks(fontsize = myfontsize)

plt.tight_layout()

plt.figure()

plt.imshow(np.abs(ch1))
plt.clim(0,0.001)

plt.show()

