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

datafolder = '/Users/boerge/skynet/molecule_computer/'
#datafolder = '/home/molecules/software/data/'


basefolder = '20191031'
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


freqs = np.genfromtxt(f_freqs, delimiter=",")
ch0 = np.genfromtxt(f_ch0, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")
ch2 = np.genfromtxt(f_ch2, delimiter=",")
ch3 = np.genfromtxt(f_ch3, delimiter=",")
ch4 = np.genfromtxt(f_ch4, delimiter=",")


# get number of averages
no_of_avg = int(len(freqs)/len(np.unique(freqs)))

print('Found ' + str(no_of_avg) + ' averages.')

freqs = av(freqs, no_of_avg)
ch0 = av(ch0, no_of_avg)
ch1 = av(ch1, no_of_avg)
ch2 = av(ch2, no_of_avg)
ch3 = av(ch3, no_of_avg)
ch4 = av(ch4, no_of_avg)



avg_freq = np.mean(freqs)

nus = 2.0*(freqs - yb_174_freq/2.0)*1e12/1e6

nus = 2*freqs

delay_in_for_loop = 60e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3




def plot_sub(k, times, nus, ch, myfontsize = 8):

    plt.subplot(2,3,k)

    plt.pcolor(times, nus, ch)

    plt.xlabel('Time (ms)', fontsize = myfontsize)
    plt.ylabel('Frequencies (MHz)', fontsize = myfontsize)
    plt.tick_params(labelsize=myfontsize,direction='in')

    plt.gca().invert_yaxis()


    plt.tight_layout()



plt.figure(figsize=(10,6))

plot_sub(1, times, nus, ch0)
plot_sub(2, times, nus, ch1)
plot_sub(3, times, nus, ch2)
plot_sub(4, times, nus, ch3)
plot_sub(5, times, nus, ch4)



plt.figure(figsize=(10,6))

plt.subplot(2,1,1)
plt.plot(times, np.sum(ch2, axis = 0))
plt.subplot(2,1,2)
plt.plot(nus, np.sum(ch2, axis = 1))

plt.show()











asd



#avg_freq = np.mean(freqs)
#avg_freq = 2*avg_freq

#nus = (freqs - yb_174_freq/2.0 )*1e12/1e6
#nus = (freqs - yb_174_freq/2.0 )*1e12/1e6




   

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


# w = k*v * cos(theta)
# v = nu * lambda/cos(theta)

velocities = -398*10**-9 * (2.0*freqs - yb_174_freq)*1e12 / np.cos(np.pi/4.0)



plt.figure(figsize=(10,6))

plt.pcolor(times, nus, ch1)
plt.colorbar()

plt.xlabel('Time (ms)', fontsize = 16)
plt.ylabel('Frequencies (MHz)', fontsize = 16)
plt.tick_params(labelsize=14,direction='in')

plt.gca().invert_yaxis()




plt.figure(figsize=(10,6))

plt.pcolor(times, nus, ch3)
plt.colorbar()

plt.xlabel('Time (ms)', fontsize = 16)
plt.ylabel('Frequencies (MHz)', fontsize = 16)
plt.tick_params(labelsize=14,direction='in')

plt.gca().invert_yaxis()


plt.xlim(0, 10)


plt.tight_layout()





plt.figure(figsize=(10,6))

plt.pcolor(times, velocities, ch3)
plt.colorbar()

plt.xlabel('Time (ms)', fontsize = 16)
plt.ylabel('Velocity (m/s)', fontsize = 16)
plt.tick_params(labelsize=14,direction='in')

#plt.gca().invert_yaxis()


plt.xlim(0, 10)
#plt.ylim(-200, 550)


plt.tight_layout()






plt.figure(figsize=(10,6))


count_min = 0.35
count_max = np.max(np.max(ch3))

#count_min = 0.75
#count_max = 3.75


ch3_plot = np.copy(ch3)
ind = np.where(ch3<count_min)

#ch3_plot[ind] = np.nan


#plt.pcolor(times, velocities, ch3_plot)
plt.contour(times, velocities, ch3_plot, np.linspace(count_min, count_max, 15))

plt.xlabel('Time (ms)', fontsize = 16)
plt.ylabel('Velocity Yb-174 (m/s)', fontsize = 16)
plt.tick_params(labelsize=14,direction='in')

#plt.gca().invert_yaxis()



vel_shift_172 = -(589.0+531.0)/2.0*1e6 * 398*10**-9 / np.cos(np.pi/4.0)
vel_shift_171 = -835.0*1e6 * 398*10**-9 / np.cos(np.pi/4.0)

yag_fire_time = 0.3

time1 = np.where(times > 1.5)[0][0] #yag_fire_time)[0][0]

plt.plot(times[time1:], 32 * 25.4e-3/(times[time1:] - yag_fire_time)/1e-3, 'r--', linewidth = 1)

plt.plot(times[time1:], vel_shift_172 + 32 * 25.4e-3/(times[time1:] - yag_fire_time)/1e-3, 'r--', linewidth = 1)

plt.plot(times[time1:], vel_shift_171 + 32 * 25.4e-3/(times[time1:] - yag_fire_time)/1e-3, 'r--', linewidth = 1)



plt.xlim(0, 10)
#plt.ylim(-200, 530)


plt.tight_layout()


plt.text(3.0, 450, 'Yb-174', fontsize = 16)
plt.text(3.0, 120, 'Yb-172/173', fontsize = 16)
plt.text(1.5, -75, 'Yb-171', fontsize = 16)



plt.savefig('ytterbium_beam.png', dpi=600)



#plt.figure(figsize=(10,6))
#
#time1 = np.where( np.abs(times - 1.75) < 0.5 )[0][0]
#time2 = np.where( np.abs(times - 4.00) < 0.5 )[0][0]
#
#vel_sum = np.mean(ch3_plot[:, time1:time2], axis = 1)
#
#plt.plot(velocities, vel_sum)
#
#plt.xlabel('Velocities (m/s)')
#plt.ylabel('Signal (a.u.)')



plt.show()

