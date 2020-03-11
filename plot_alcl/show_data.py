import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
from fit_yb import *
from fit_mo import *


c = 299792458
yb_174_freq = 751.52653349 # in THz


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

#datafolder = '/Users/boerge/skynet/molecule_computer/'
datafolder = '/home/molecules/software/data/'


basefolder = '20200309'
#time_stamp = '123451'
time_stamp = sys.argv[1]

#basefolder = '20191004'
#time_stamp = '171345'



basefilename = datafolder + basefolder + '/' + basefolder + '_'


f_freqs = basefilename + time_stamp + '_set_points'
f_ch1 = basefilename + time_stamp + '_ch0_arr'
f_ch2 = basefilename + time_stamp + '_ch1_arr'
f_ch3 = basefilename + time_stamp + '_ch2_arr'
f_ch3 = basefilename + time_stamp + '_ch3_arr'


freqs = np.genfromtxt(f_freqs, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")
ch2 = np.genfromtxt(f_ch2, delimiter=",")
ch3 = np.genfromtxt(f_ch3, delimiter=",")


# get number of averages
no_of_avg = int(len(freqs)/len(np.unique(freqs)))

print('Found ' + str(no_of_avg) + ' averages.')

freqs = av(freqs, no_of_avg)
ch1 = av(ch1, no_of_avg)
ch2 = av(ch2, no_of_avg)
ch3 = av(ch3, no_of_avg)


#avg_freq = np.mean(freqs)
#avg_freq = 2*avg_freq

#nus = (freqs - yb_174_freq/2.0 )*1e12/1e6
#nus = (freqs - yb_174_freq/2.0 )*1e12/1e6

avg_freq = np.mean(freqs)

nus = 2.0*(freqs - yb_174_freq/2.0)*1e12/1e6




delay_in_for_loop = 60e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3

   

cut_time1 = 10.0
cut_time2 = 12.0

ch1_start = np.where( np.abs(times - cut_time1) < 0.5 )[0][0]
ch1_end = np.where( np.abs(times - cut_time2) < 0.5 )[0][0]

# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
    ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
    ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])
    ch3[k, :] = ch3[k, :] - np.mean(ch3[k, -offset_avg_points:-1])


ind = np.where(ch1 > 0.005)

ch1[ind] = 0.0


plt.figure()
plt.subplot(2,1,1)
plt.plot(np.mean(ch1[:, :], axis = 0))

c1 = 189
c2 = 200

plt.subplot(2,1,2)
plt.plot(np.mean(ch1[:, c1:c2], axis = 1))
plt.plot(np.mean(ch1[:, 10:20], axis = 1))

plt.figure()
plt.pcolor(times, nus, ch1)
plt.colorbar()

plt.show()

