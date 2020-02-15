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


#basefolder = str(my_today.strftime('%Y%m%d')) # 20190618
basefolder = '20190628'

#basefolder = '20190910'
#basefolder = '20190627'

#basefilename = datafolder + basefolder + '/' + basefolder + '_' # 20190618_105557
basefilename = datafolder + basefolder + '/' + basefolder + '_'

if len(sys.argv)>1:
    time_stamp = sys.argv[1]
else:
    # get latest time stamp
    all_files = np.sort(glob.glob(basefilename + "*"))
    #print(all_files)
    time_stamp = all_files[-1].split('_')[1]


# molybdenum data
time_stamp = '144843'


f_freqs = basefilename + time_stamp + '_freqs'
f_ch1 = basefilename + time_stamp + '_ch1'
f_ch2 = basefilename + time_stamp + '_ch2'
f_ch3 = basefilename + time_stamp + '_ch3'


freqs = np.genfromtxt(f_freqs, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")
ch2 = np.genfromtxt(f_ch2, delimiter=",")

#ch3 = np.genfromtxt(f_ch3, delimiter=",")
ch3 = ch2


# get number of averages
no_of_avg = int(len(freqs)/len(np.unique(freqs)))

print('Found ' + str(no_of_avg) + ' averages.')

freqs = av(freqs, no_of_avg)
ch1 = av(ch1, no_of_avg)
ch2 = av(ch2, no_of_avg)
ch3 = av(ch3, no_of_avg)


avg_freq = np.mean(freqs)
avg_freq = 2*avg_freq

#nus = (freqs - yb_174_freq/2.0 )*1e12/1e6
#nus = (freqs - yb_174_freq/2.0 )*1e12/1e6

avg_freq = np.mean(freqs)

nus = 2*(freqs - avg_freq)*1e12/1e6





delay_in_for_loop = 60e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3


# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
        ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
        ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])
        ch3[k, :] = ch3[k, :] - np.mean(ch3[k, -offset_avg_points:-1])

    
   


#freq1_ind = np.where( np.abs(nus - -500.0) < 30.0 )[0][0]
#freq2_ind = np.where( np.abs(nus - 0.0) < 30.0 )[0][0]
#
#
#plt.figure()
#plt.plot(times, np.mean(ch3[freq1_ind:freq2_ind, :], axis = 0))
#plt.xlabel('Time (ms)')
#plt.ylabel('Signal (a.u.)')
#
#plt.figure()
#plt.subplot(2,1,1)
#plt.pcolor(times, nus, ch3)
#plt.xlabel('Time (ms)')
#plt.ylabel('Frequency (MHz) + ' + str(avg_freq) + ' THz')
#
#plt.axhline(nus[freq1_ind], color = 'r')
#plt.axhline(nus[freq2_ind], color = 'r')
#
#plt.subplot(2,1,2)
#plt.pcolor(times, nus, ch1)
#plt.xlabel('Time (ms)')
#plt.ylabel('Frequency (MHz) + ' + str(avg_freq) + ' THz')
#
#plt.axhline(nus[freq1_ind], color = 'r')
#plt.axhline(nus[freq2_ind], color = 'r')






cut_time1 = 0.2
cut_time2 = 3.0

ch1_start = np.where( np.abs(times - cut_time1) < 0.5 )[0][0]
ch1_end = np.where( np.abs(times - cut_time2) < 0.5 )[0][0]

# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
    ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
    ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])
    ch3[k, :] = ch3[k, :] - np.mean(ch3[k, -offset_avg_points:-1])





spectrum = np.mean(ch1[:, ch1_start:ch1_end], axis = 1)


x_fit = 0
y_fit = 0
(x_fit, y_fit, result) = fit_mo(nus, spectrum)




# shifting the zero point in the plot to Mo-96
my_shift = result.params['x_offset'].value # in MHz

nus = nus + my_shift
x_fit = x_fit + my_shift
avg_freq = avg_freq - my_shift*1e6/1e12



for k in result.params.keys():
    print(str(k) + ' = ' + str(result.params[k].value))




fig = plt.figure(figsize=(10,6))

moly   = np.array([100, 98, 97, 96, 95, 94]) # isotopes
vshift = np.array([-0.5,-0.8,-0.7,-0.7,-0.6,-0.6])*-1
hshift = np.array([20,-160,-130,10,-140,20])


# moly isotope shifts
text_freqs = my_shift + np.array([-0.7945,-0.2938,-0.0180,0,+0.3028,+0.4107,0.9144]) * 1000 # in MHz


for k in range(len(moly)):
    plt.axvline(text_freqs[k] - result.params['x_offset'], linestyle =  '--',linewidth=1.6,label='Mo'+str(moly[k]))
    plt.text(text_freqs[k] - result.params['x_offset'] + hshift[k], vshift[k]   , 'Mo ' + str(moly[k]),fontsize=16)


# plot data and fit

plt.plot(x_fit, -1*y_fit, 'k-', linewidth = 3)
plt.scatter(nus, -1*spectrum, color = 'r', edgecolor = 'k', s = 100)
plt.plot(nus, -1*spectrum, color = 'r', linestyle = '-') 


plt.xlabel('Frequency (MHz) + ' + "{0:2.6f}".format(avg_freq) + ' THz',fontsize=16)
plt.ylabel('Signal (a.u)',fontsize=16)
plt.tick_params(labelsize=14,direction='in')
#plt.xlim(-760,760)
plt.xlim(np.min(nus), np.max(nus))


plt.tight_layout()
#plt.legend(fontsize=14)


plt.savefig('molybdenum.png', dpi=600)

plt.show()


