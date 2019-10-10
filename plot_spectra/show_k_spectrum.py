import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
from fit_k import *


c = 299792458
#K_39_freq = 389.286058716 # D1-line in THz
K_39_freq = 391.01617003 - 173.1e-6 + 14.4e-6#  D2-line in THz, F=2 -> F'=3

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

datafolder = '/Users/boerge/skynet/K_Tests/'

datafolder = '/Users/boerge/software/Prehistoric-Data-Acquisition/'


# molybdenum data
time_stamp = sys.argv[1]



#basefolder = str(my_today.strftime('%Y%m%d')) # 20190618
basefolder = '20190620'

basefilename = datafolder + '2019-05-02-'

f_freqs = basefilename + time_stamp + '_s' + '.csv'
f_ch1 = basefilename + time_stamp + '_1' + '.csv'
f_ch2 = basefilename + time_stamp + '_2' + '.csv'

basefilename = datafolder

f_freqs = basefilename + 'setpts.txt'
f_ch1 = basefilename + 'ch1.txt'
f_ch2 = basefilename + 'ch2.txt'


freqs = np.genfromtxt(f_freqs, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")
ch2 = np.genfromtxt(f_ch2, delimiter=",")



avg_freq = K_39_freq
nus = (freqs - avg_freq)*1e12/1e6





delay_in_for_loop = 60e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3





cut_time1 = 50.0
cut_time2 = 55.0

ch1_start = np.where( np.abs(times - cut_time1) < 0.5 )[0][0]
ch1_end = np.where( np.abs(times - cut_time2) < 0.5 )[0][0]

# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
    ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
    ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])



## sort data if random scan was used
#ind = np.argsort(freqs)
#
#freqs = freqs[ind]
#ch1 = ch1[ind, :]
#ch2 = ch2[ind, :]




plt.figure()
plt.subplot(2,1,1)
plt.pcolor(ch1)
plt.subplot(2,1,2)
plt.pcolor(ch2)

spectrum = np.mean(ch1[:, ch1_start:ch1_end], axis = 1)




x_fit = 0
y_fit = 0
(x_fit, y_fit, result) = fit_k(nus, spectrum)




# shifting the zero point in the plot to K-39
my_shift = result.params['x_offset'].value # in MHz

nus = nus + my_shift
x_fit = x_fit + my_shift
avg_freq = avg_freq - my_shift*1e6/1e12



for k in result.params.keys():
    print(str(k) + ' = ' + str(result.params[k].value))




fig = plt.figure(figsize=(10,6))



pot    = np.array([39, 39]) # isotopes
vshift = np.array([22, 22, 25, 25]) # vertical shift of the text
hshift = np.array([ +50,  +30, 0,     0]) # horizontal shift of the text




# pottasium isotope shifts
text_freqs  = my_shift + +173.1-14.4+np.array([+288.6-16.1, -173.1+14.4, +236.2+158.8-8.4, +236.2-95.3+8.4])


for k in range(len(pot)):
    plt.axvline(text_freqs[k] - result.params['x_offset'], linestyle =  '--',linewidth=1.6,label='K'+str(pot[k]))
    plt.text(text_freqs[k] - result.params['x_offset'] + hshift[k], vshift[k]   , 'K ' + str(pot[k]),fontsize=16)


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


plt.savefig('pottassium.png', dpi=600)

plt.show()

