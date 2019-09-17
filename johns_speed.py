import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
# from fit_yb import *
# from fit_mo import *

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
#datafolder = '/Users/boerge/Github/Data/'
datafolder = '/home/molecules/software/data/'
basefolder = str(my_today.strftime('%Y%m%d')) # 20190618

basefolder = '20190806'

basefilename = datafolder + basefolder + '/' + basefolder + '_' # 20190618_105557


if len(sys.argv)>1:
    time_stamp = sys.argv[1]
else:
    # get latest time stamp
    all_files = np.sort(glob.glob(basefilename + "*"))
    time_stamp = all_files[-1].split('_')[1]

f_freqs = basefilename + time_stamp + '_freqs'
f_ch1 = basefilename + time_stamp + '_ch1'
f_ch2 = basefilename + time_stamp + '_ch2'
f_ch3 = basefilename + time_stamp + '_ch3'



cut_time1 = 0 # ms
cut_time2 = 2 # ms
freq_cut = 50 # MHz

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

avg_freq = np.mean(freqs)

nus = (freqs - avg_freq)*1e12/1e6

delay_in_for_loop = 50e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (8.97e-6 + delay_in_for_loop) / 1e-3

ch1_start = np.where( np.abs(times - cut_time1) < 0.5 )[0][0]
ch1_end = np.where( np.abs(times - cut_time2) < 0.5 )[0][0]
freq_ind = np.where( np.abs(nus - freq_cut) < 30.0 )[0][0]


# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
    ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
    ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -offset_avg_points:-1])

# scaling the frequency axis to the blue
nus = 2*nus
avg_freq = 2*avg_freq


imin = 25
imax = 45

f0 = 375.76354
adj_f = 2*freqs




spds = 0


plt.figure()
plt.imshow(ch3)
plt.show()

# spectrum = np.mean(ch1[:, ch1_start:ch1_end], axis = 1)

# (x_fit, y_fit, result) = fit_yb(nus, spectrum)




# # plotting

# plt.figure()

# # ch1

# plt.subplot(3,2,1)
# plt.pcolor(nus, times, np.transpose(ch1))#, aspect = 'auto')
# plt.xlabel('Frequency (MHz) + ' + str(avg_freq) + ' THz')
# plt.ylabel('Time (ms)')

# plt.title('Scan #: ' + time_stamp)

# plt.axhline(times[ch1_start], linewidth = 1, color = 'r', linestyle = '--')
# plt.axhline(times[ch1_end], linewidth = 1, color = 'r', linestyle = '--')
# plt.axvline(nus[freq_ind], linewidth = 1, color = 'r', linestyle = '--')

# plt.subplot(3,2,3)
# plt.plot(nus, spectrum)
# plt.plot(x_fit, y_fit)
# plt.xlabel('Frequency (MHz) + ' + str(avg_freq) + ' THz')

# plt.subplot(3,2,5)
# plt.plot(times, ch1[freq_ind, :])
# plt.xlabel('Time (ms)')


# # ch2
# plt.subplot(3,2,2)
# plt.pcolor(nus, times, np.transpose(ch2))#, aspect = 'auto')
# plt.xlabel('Frequency (MHz) + ' + str(avg_freq) + ' THz')
# plt.ylabel('Time (ms)')

# ch2_start = 0
# ch2_end = 1

# plt.axhline(times[ch2_start], linewidth = 1, color = 'r', linestyle = '--')
# plt.axhline(times[ch2_end], linewidth = 1, color = 'r', linestyle = '--')
# plt.axvline(nus[freq_ind], linewidth = 1, color = 'r', linestyle = '--')

# plt.subplot(3,2,4)
# plt.plot(nus, np.mean(ch2[:, ch2_start:ch2_end], axis = 1))
# plt.xlabel('Frequency (MHz) + ' + str(avg_freq) + ' THz')

# plt.subplot(3,2,6)
# plt.plot(times, ch2[freq_ind, :])
# plt.xlabel('Time (ms)')



# plt.tight_layout()





# # shifting the zero point in the plot to Yb-174
# yb_174_freq = 751.52653349 # in THz

# my_shift = result.params['x_offset'].value # in MHz

# nus = nus + my_shift
# x_fit = x_fit + my_shift
# avg_freq = avg_freq - my_shift*1e6/1e12

# print('Wavemeter offset ' + str((yb_174_freq - avg_freq) * 1e6) + ' MHz')


# plt.figure()

# plt.scatter(nus, -1*spectrum, color = 'r', edgecolor = 'b')
# plt.plot(nus, -1*spectrum, color = 'm', linestyle = '-') 
# plt.plot(x_fit, -1*y_fit, 'k-')
# plt.xlabel('Measured Frequency (MHz) + ' + "{0:2.6f}".format(avg_freq) + ' THz',fontsize=16)
# plt.ylabel('Signal (a.u)',fontsize=16)
# plt.tick_params(labelsize=14,direction='in')
# plt.xlim(np.min(nus), np.max(nus))

# freqs  = my_shift + np.array([-508.89, 0 , -250.78, 589.75, 531.11, 835.19, 1153.68, 1190.36, 1888.80])
# yb     = np.array([176    ,174,     173,    173,    172,    171,     171,     170,     168]) # isotopes
# vshift = np.array([2.7    ,0.1,     2.0,    0.6,   0.35,   0.25,     2.0,     2.0,       0]) # vertical shift of the text
# hshift = np.array([ 15  , -130,    -140,     15,   -140,     20,    -140,      15,       0]) # horizontal shift of the text
# for k in result.params.keys():
#     print(str(k) + ' = ' + str(result.params[k].value))

# for k in range(len(yb)):
#     plt.axvline(freqs[k] - result.params['x_offset'], linestyle =  '--',linewidth=1.6,label='Yb'+str(yb[k]))
#     plt.text(freqs[k] - result.params['x_offset'] + hshift[k], vshift[k]   , 'Yb ' + str(yb[k]),fontsize=16)



# plt.show()


