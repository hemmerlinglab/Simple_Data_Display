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
#basefolder = '20200121'
basefolder = '20200211'
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


plt.figure()
plt.pcolor(times, freqs, ch1)




spectrum = np.mean(ch1[:, ch1_start:ch1_end], axis = 1)


x_fit = 0
y_fit = 0
(x_fit, y_fit, result) = fit_yb(nus, spectrum)




# shifting the zero point in the plot to Mo-96
my_shift = result.params['x_offset'].value # in MHz

nus = nus + my_shift
x_fit = x_fit + my_shift
avg_freq = avg_freq - my_shift*1e6/1e12



for k in result.params.keys():
    print(str(k) + ' = ' + str(result.params[k].value))




fig = plt.figure(figsize=(10,6))



yb     = np.array([176    , 174,  173,    173,    172,    171,     171,     170]) # isotopes
vshift = np.array([1.35   ,1.35,  1.2,    1.2,   1.35,   1.25,    1.35,    1.35]) # vertical shift of the text
hshift = np.array([ 15  ,    50, -190,     15,   -200,     10,    -190,      15]) # horizontal shift of the text




# moly isotope shifts
text_freqs  = my_shift + np.array([-508.89, 0 , -250.78, 589.75, 531.11, 835.19, 1153.68, 1190.36, 1888.80])


for k in range(len(yb)):
    plt.axvline(text_freqs[k] - result.params['x_offset'], linestyle =  '--',linewidth=1.6,label='Yb'+str(yb[k]))
    plt.text(text_freqs[k] - result.params['x_offset'] + hshift[k], vshift[k]   , 'Yb ' + str(yb[k]),fontsize=16)


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


plt.savefig('ytterbium.png', dpi=600)




plt.show()

