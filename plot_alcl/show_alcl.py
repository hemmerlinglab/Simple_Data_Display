import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
from fit_line import *

c = 299792458


# hlp is helper variable
def plot_sub(k, times, nus, ch, myfontsize = 8):

    plt.subplot(3,2,k)

    plt.pcolor(times, nus, ch)

    plt.xlabel('Time (ms)', fontsize = myfontsize)
    plt.ylabel('Frequencies (MHz)', fontsize = myfontsize)
    plt.tick_params(labelsize=myfontsize,direction='in')

    plt.gca().invert_yaxis()


    plt.tight_layout()

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
#datafolder = '/home/molecules/software/data/'
datafolder = '/Users/boerge/software/data/Data/molecule_computer/'


basefolder = '20200211'
time_stamp = sys.argv[1]




basefilename = datafolder + basefolder + '/' + basefolder + '_'


f_freqs = basefilename + time_stamp + '_set_points'
f_ch0 = basefilename + time_stamp + '_ch0_arr'
f_ch1 = basefilename + time_stamp + '_ch1_arr'
f_ch2 = basefilename + time_stamp + '_ch2_arr'
f_ch3 = basefilename + time_stamp + '_ch3_arr'
f_ch4 = basefilename + time_stamp + '_ch4_arr'

f_ch0s = basefilename + time_stamp + '_ch0_arr_slow'
f_ch1s = basefilename + time_stamp + '_ch1_arr_slow'
f_ch2s = basefilename + time_stamp + '_ch2_arr_slow'
f_ch3s = basefilename + time_stamp + '_ch3_arr_slow'
f_ch4s = basefilename + time_stamp + '_ch4_arr_slow'



freqs = np.genfromtxt(f_freqs, delimiter=",")
ch0 = np.genfromtxt(f_ch0, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")
ch2 = np.genfromtxt(f_ch2, delimiter=",")
ch3 = np.genfromtxt(f_ch3, delimiter=",")
ch4 = np.genfromtxt(f_ch4, delimiter=",")

ch0s = np.genfromtxt(f_ch0s, delimiter=",")
ch1s = np.genfromtxt(f_ch1s, delimiter=",")
ch2s = np.genfromtxt(f_ch2s, delimiter=",")
ch3s = np.genfromtxt(f_ch3s, delimiter=",")
ch4s = np.genfromtxt(f_ch4s, delimiter=",")



# get number of averages
no_of_avg = int(len(freqs)/len(np.unique(freqs)))

print('Found ' + str(no_of_avg) + ' averages.')

freqs = av(freqs, no_of_avg)
ch0 = av(ch0, no_of_avg)
ch1 = av(ch1, no_of_avg)
ch2 = av(ch2, no_of_avg)
ch3 = av(ch3, no_of_avg)
ch4 = av(ch4, no_of_avg)

ch0s = av(ch0s, no_of_avg)
ch1s = av(ch1s, no_of_avg)
ch2s = av(ch2s, no_of_avg)
ch3s = av(ch3s, no_of_avg)
ch4s = av(ch4s, no_of_avg)




offset_freq = 382.0959

#nus = 2*freqs

delay_in_for_loop = 100e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3




time_cut1 = 98
time_cut2 = 98+20


color_plot = ch0[:, time_cut1:time_cut2]

color_plot -= np.min(np.min(color_plot))

color_plot /= np.max(np.max(color_plot))


plt.figure()
plt.pcolor(times[time_cut1:time_cut2], 3*freqs, color_plot)
plt.xlabel('Time (us)')
plt.ylabel('Relative UV frequency (MHz)')

plt.colorbar()

plt.figure()

signal = np.mean(ch0[:, time_cut1:time_cut2], axis = 1)
signal = signal/np.max(np.abs(signal))

(xfit, yfit, result) = fit_line(freqs, signal)

cnt_peak = result.params['x0'].value

plt.plot(3*freqs, 1 - signal)
plt.plot(3*xfit, 1 - yfit, 'r-')

plt.ylim(0, 0.5)

plt.ylabel('Absorption signal (a.u.)')
plt.xlabel('Frequency (MHz)')

abs_cnt_peak = 3 * ((offset_freq * 1e12 + cnt_peak * 1e6)/1e12) # in THz

abs_cnt_peak_wavenumber = abs_cnt_peak * 1e12/100.0/c

plt.text(-400, 0.4, "Center peak frequency: \n\n     {0:9.6f} THz \n = {1:9.4f} 1/cm".format(abs_cnt_peak, abs_cnt_peak_wavenumber))



plt.show()

