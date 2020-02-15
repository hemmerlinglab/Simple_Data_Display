import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
from fit_k import *


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


datafolder = '/home/molecules/software/data/'


# molybdenum data


basefolder = '20200213'

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

#ch1 = subtract_noise(ch1,ch4)


nus = freqs*1.5


delay_in_for_loop = 100e-6
#delay_in_for_loop = 50e-6

no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3





cut_time1 = 10.0
cut_time2 = 11.0

ch1_start = np.where( np.abs(times - cut_time1) < 0.5 )[0][0]
ch1_end = np.where( np.abs(times - cut_time2) < 0.5 )[0][0]

# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
    ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
    ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])
    ch4[k, :] = ch4[k, :] - np.mean(ch4[k, -offset_avg_points:-1])




plt.figure()
#plt.pcolor(data1[:, 98:300])
plt.pcolor(times[98:300],nus,ch1[:, 98:300])

#plt.pcolor(times[182:300],nus,ch1[:, 182:300])

lines = {
124050:382.08140,
130453:382.08615,
132440:382.09095,
134013:382.09340,
135425:382.09575,
140938:382.09810
}

plt.xlabel('Time (ms)')
plt.ylabel('Frequency - {} (THz)'.format(lines[int(time_stamp)]*3))

plt.title('{} THz ({})'.format(lines[int(time_stamp)]*3,time_stamp))

#plt.figure()

#plt.plot(np.mean(ch1[:, 99:120], axis = 1))
#plt.plot(np.mean(ch1[:, 182:200], axis = 1))
#plt.plot(np.mean(data1[:, 99:120], axis = 1))

#plt.figure()

#plt.plot(np.mean(ch1[:, 90:150], axis = 0))
#plt.plot(np.mean(ch1[:, 170:200], axis = 0))
#plt.ylim([-0.005, 0.005])


plt.show()

