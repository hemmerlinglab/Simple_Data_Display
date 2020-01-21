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

#basefolder = str(my_today.strftime('%Y%m%d')) # 20190618
basefolder = '20200117'

basefilename = datafolder + basefolder + '/' + basefolder + '_'

if len(sys.argv)>1:
    time_stamp = sys.argv[1]
else:
    # get latest time stamp
    all_files = np.sort(glob.glob(basefilename + "*"))
    #print(all_files)
    time_stamp = all_files[-1].split('_')[1]


# molybdenum data
#time_stamp = '122803'


f_freqs = basefilename + time_stamp + '_set_points'
f_ch1 = basefilename + time_stamp + '_ch0_arr'
f_ch2 = basefilename + time_stamp + '_ch1_arr'


freqs = np.genfromtxt(f_freqs, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")
ch2 = np.genfromtxt(f_ch2, delimiter=",")


# get number of averages
no_of_avg = int(len(freqs)/len(np.unique(freqs)))

print('Found ' + str(no_of_avg) + ' averages.')

freqs = av(freqs, no_of_avg)
ch1 = av(ch1, no_of_avg)
ch2 = av(ch2, no_of_avg)


avg_freq = np.mean(freqs)
avg_freq = 2*avg_freq

#nus = (freqs - yb_174_freq/2.0 )*1e12/1e6
#nus = (freqs - yb_174_freq/2.0 )*1e12/1e6

avg_freq = np.mean(freqs)

nus = freqs - yb_174_freq*1e12/1e6*0


delay_in_for_loop = 100e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3


# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
        ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
        ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])

    
   

cut_time1 = 10.0
cut_time2 = 12.0

ch1_start = np.where( np.abs(times - cut_time1) < 0.5 )[0][0]
ch1_end = np.where( np.abs(times - cut_time2) < 0.5 )[0][0]

# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
    ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
    ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])





spectrum = np.mean(ch1[:, ch1_start:ch1_end], axis = 1)


plt.figure()
#plt.plot(spectrum)
plt.subplot(2,2,1)
plt.pcolor(times, freqs, ch1)
plt.subplot(2,2,3)
plt.pcolor(times, freqs, ch2)

plt.subplot(2,2,2)
plt.plot(times, np.mean(ch1, axis = 0))


plt.figure()

plt.plot(ch1)

plt.show()

