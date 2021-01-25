import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
from helper import *
from fit_rb import *
from scipy.signal import find_peaks
from fit_lorentz import fit_lorentz


my_date = '20200616'
my_time = '110953'

my_date = '20200618'
my_time = '104630'




#folder = '/home/molecules/software/data/'
folder = '/Users/boerge/Software/offline_data/'

(freqs, signal, wavemeter_freqs) = get_data(my_date, my_time, datafolder = folder)

cut1 = 0
cut2 = 180

#cut1 = 0
#cut2 = len(freqs)

freqs = freqs[cut1:cut2]
#wavemeter_freqs = wavemeter_freqs[cut1:cut2]

n_signal = np.zeros([3, cut2-cut1])

n_signal[0, :] = signal[0, cut1:cut2]
n_signal[1, :] = signal[1, cut1:cut2]
n_signal[2, :] = signal[2, cut1:cut2]

signal = n_signal



# get Rubidium frequencies
rb_transitions = calculate_Rb_transitions()



cnt_freq = 384.230484468e12 # S1/2 - P3/2
myfontsize = 16

x_plot = (freqs - cnt_freq)/1e9


guess_wavemeter_offset = -10.0e6


# calculate difference signal
diff = signal[0] - signal[1]
diff -= np.min(diff)




plt.figure()
plt.subplot(2,1,1) 
plt.title('Scan ' + my_date + "_" + my_time)


plt.plot(x_plot, diff, '.-')


plt.xlim(np.min(x_plot), np.max(x_plot))

line_centers_fit = []
line_centers_true = []
for k in range(len(rb_transitions)):

    peak_ind = np.where( np.abs(x_plot - (rb_transitions[k][0] - cnt_freq + guess_wavemeter_offset)/1e9) < 0.005 )[0][0]
    
    (xf, yf, peak_result, res) = fit_lorentz(x_plot, diff, peak_ind, dn = 5)

    plt.plot(xf, yf, 'k')

    line_centers_fit.append(peak_result.params['cnt0'].value*1e9 + cnt_freq)
    
    # has to have some number of peaks and was given
    line_centers_true.append(rb_transitions[k][0])

    plt.axvline((rb_transitions[k][0] - cnt_freq)/1e9, ls = '--', color = 'k', label = 'True')
    
    plt.text((rb_transitions[k][0] - cnt_freq)/1e9 + 0.01, 0.0, rb_transitions[k][1], rotation = 90, fontsize = 8)
    
    plt.axvline((line_centers_fit[-1] - cnt_freq)/1e9, ls = '--', color = 'r', label = 'Fit')

line_centers_fit = np.sort(np.array(line_centers_fit))
line_centers_true = np.sort(np.array(line_centers_true))

#plt.legend()

plt.xlim(-2.7, -2.2)

# plot true value vs fitted line centers

plt.subplot(2,1,2) 
plt.plot((line_centers_true - cnt_freq)/1e9, (line_centers_fit - line_centers_true)/1e6, 'o-')

plt.xlabel('Frequency Offset (GHz) + {0:2.6f} THZ'.format(cnt_freq/1e12));
plt.ylabel('Fit - True Frequency (MHz)');

plt.xlim(np.min(x_plot), np.max(x_plot))

plt.xlim(-2.7, -2.2)

# fitted average offset
mean_wavemeter_offset = np.mean(line_centers_fit - line_centers_true)/1e6

plt.axhline(mean_wavemeter_offset, ls = '--')
plt.text(-2.5, 0.05 + mean_wavemeter_offset, "{0:2.2f} MHz".format(mean_wavemeter_offset))


plt.savefig(my_date + '_' + my_time + '_Rb_calibration.png', format = 'png')

plt.show()




