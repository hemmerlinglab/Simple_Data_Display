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
datafolder = '/home/molecules/software/data/'


basefolder = '20191107'
#time_stamp = '123451'
time_stamp = sys.argv[1]

#basefolder = '20191004'
#time_stamp = '171345'



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


plt.figure()
plt.plot(np.sum(ch0[:, 180:185], axis = 1))
plt.show()
asd

avg_freq = np.mean(freqs)

nus = 2.0*(freqs - yb_174_freq/2.0)*1e12/1e6

#nus = 2*freqs

delay_in_for_loop = 100e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3



# change nus to vel

#vel = -(2.0*freqs*1e6 + 2.0*375.763216e12 - yb_174_freq*1e12) * 398e-9 / np.cos(np.pi/4.0)
vel = -(freqs*2.0 - yb_174_freq)*1e12 * 398e-9 #/ np.cos(np.pi/4.0)


plt.figure(figsize=(10,6))
plt.pcolor(times, vel, ch0)

plt.figure(figsize=(10,6))
plt.pcolor(times, vel, ch2)

dist = 31*25.4e-3

yag_fire = 9.567e-3

isotope_shifts = np.array([-508.89, 0, -250.78, (531.11 + 589.75)/2.0, 835.19, (1153.68 + 1190.36)/2.0, 1888.80])

vel_shifts = -(isotope_shifts)*1e6 * 398e-9 #/ np.cos(np.pi/4.0)

for k in range(len(vel_shifts)):
    if k == 1:
        plt.plot(times, vel_shifts[k] + dist/(times*1e-3 - yag_fire), 'k--')
    else:
        plt.plot(times, vel_shifts[k] + dist/(times*1e-3 - yag_fire), 'r--')

plt.ylim(np.min(vel), np.max(vel))
plt.xlim(yag_fire/1e-3, 20)
plt.show()

#asd

plt.figure(figsize=(10,6))

plot_sub(1, times, nus, ch0)
plot_sub(3, times, nus, ch0s)

cmax = 1.5

plot_sub(2, times, nus, ch2)
plt.colorbar()
#plt.clim(0,cmax)

plot_sub(4, times, nus, ch2s)
plt.colorbar()
#plt.clim(0,cmax)

#plot_sub(2, times, nus, ch1s)
#plot_sub(4, times, nus, ch3s)
#plot_sub(5, times, nus, ch4s)

plot_sub(5, times, nus, ch0s - ch0)
plot_sub(6, times, nus, ch2s - ch2)
plt.colorbar()
plt.clim(-0.25*cmax, 0.25*cmax)

plt.figure(figsize=(10,6))

plt.subplot(2,1,1)
plt.plot(times, np.sum(ch2s - ch2, axis = 0))
#plt.subplot(2,1,2)
#plt.plot(nus, np.sum(ch2, axis = 1))


plt.show()


