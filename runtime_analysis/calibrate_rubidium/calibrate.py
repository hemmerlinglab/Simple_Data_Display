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



my_date = '20200713'
#my_time = '145145'

my_time = sys.argv[1]

#my_date = '20200619'
#my_time = '095302'

#my_date = '20200618'
#my_time = '104630'

folder = '/home/molecules/software/data/'
#folder = '/Users/boerge/Software/offline_data/'

(freqs, signal, wavemeter_freqs) = get_data(my_date, my_time, datafolder = folder)

#cut1 = 0
#cut2 = 90

#cut1 = 90
#cut2 = 160

cut1 = 0
cut2 = 180

cut1 = 0
cut2 = len(freqs)

freqs = freqs[cut1:cut2]
wavemeter_freqs = wavemeter_freqs[cut1:cut2]

n_signal = np.zeros([3, cut2-cut1])

n_signal[0, :] = signal[0, cut1:cut2]
n_signal[1, :] = signal[1, cut1:cut2]
n_signal[2, :] = signal[2, cut1:cut2]

signal = n_signal


rb_transitions = calculate_Rb_transitions()

# fit spectrum

(xfit, yfit, fit_result, residuals) = fit_rb(freqs, signal[0], rb_transitions)

wavemeter_offset = fit_result.params['x_offset'].value
print('Wavemeter offset: {0:2.6} MHz'.format(wavemeter_offset/1e6))




cnt_freq = 384.230484468e12 # S1/2 - P3/2
cnt_freq = 377.107e12 # S1/2 - P1/2
myfontsize = 16

x_plot = (freqs - cnt_freq)/1e9

plt.figure(figsize = (20,6))
plt.subplot(2,1,1)

plt.plot(x_plot, signal[0], 'r.-')
plt.plot(x_plot, signal[1], 'g-')

plt.plot((xfit - cnt_freq)/1e9, yfit, 'k-')

plt.xlabel("Frequency (GHz) - {0:2.6f} THz".format(cnt_freq/1e12), fontsize = myfontsize)
plt.ylabel('Absorption Signal (a.u)', fontsize = myfontsize)


plt.tick_params(labelsize=14, direction='in')



dv = 1.0/len(rb_transitions)
for k in range(len(rb_transitions)):

    plt.axvline((rb_transitions[k][0] + wavemeter_offset - cnt_freq)/1e9, ls = rb_transitions[k][2], color = 'r')
    #plt.axvline((rb_transitions[k][0] - cnt_freq)/1e9, ls = rb_transitions[k][2], color = 'k')

    plt.text((rb_transitions[k][0] + wavemeter_offset - cnt_freq)/1e9 + 0.01, 1.0, rb_transitions[k][1], rotation = 90, fontsize = 8)
    


plt.text(-2.2, np.min(yfit), 'Wavemeter offset: {0:2.6} MHz'.format(wavemeter_offset/1e6), fontsize = myfontsize)

plt.xlim(min(x_plot), max(x_plot))

plt.title('Rb Calibration Scan - ' + my_date + '_' + my_time)

plt.tight_layout()


plt.subplot(2,1,2)
plt.plot(x_plot, (wavemeter_freqs - freqs)/1e6)
plt.xlabel("Frequency (GHz) - {0:2.6f} THz".format(cnt_freq/1e12), fontsize = myfontsize)
plt.ylabel('Set - Act freq (MHz)', fontsize = myfontsize)
plt.tick_params(labelsize=14, direction='in')
plt.xlim(min(x_plot), max(x_plot))

plt.tight_layout()

plt.savefig(my_date + '_' + my_time + '_Rb_calibration.png', format = 'png')




##########################
## Subtract background
##########################
#
#plt.figure()
#
#y_offset = fit_result.params['y_offset']
#a0 = fit_result.params['ab0']
#a1 = fit_result.params['ab1']
#
#cnt0 = fit_result.params['cnt0']
#cnt1 = fit_result.params['cnt1']
#w0 = fit_result.params['w0']
#w1 = fit_result.params['w1']
#x_offset = fit_result.params['x_offset']
#
#x_fg = freqs
#
#x_fg = np.linspace(min(freqs)-1e9,max(freqs)+1e9,1000)
#
#gauss_f  = y_offset 
#gauss_f += -a0 * np.exp( -(x_fg - cnt0 - x_offset)**2/(2.0*w0**2) )
#gauss_f += -a1 * np.exp( -(x_fg - cnt1 - x_offset)**2/(2.0*w1**2) )
# 
##y_fg = signal[0] - gauss_f 
#
#
##x_plot = x_fg
#
#
#x_fg = (x_fg - cnt_freq)/1e9
#
#plt.plot(x_plot, signal[0])
#plt.plot(x_plot, signal[1])
#plt.plot(x_fg, gauss_f)





if True:


    #########################
    # Plot residuals
    #########################
    #
    #plt.figure()
    #
    #plt.plot((freqs - cnt_freq)/1e9, residuals, 'k-')
    #plt.xlabel("Frequency (GHz) - {0:2.6f} THz".format(cnt_freq/1e12), fontsize = myfontsize)
    #plt.xlim(min(x_plot), max(x_plot))
    #plt.ylabel('Residuals')
    #
    #plt.tight_layout()
    #
    #
    #
    plt.figure()
    plt.subplot(2,1,1) 
    diff = signal[0] - signal[1]

    diff -= np.min(diff)

    #diff = signal[0]

    #peaks, _ = find_peaks(diff, width = 5, distance = 5)
    #peaks, _ = find_peaks(diff, distance = 15)
    
    plt.plot(x_plot, diff, '.-')
    
    #plt.plot(x_plot[peaks], diff[peaks], "x")
    
    plt.xlim(np.min(x_plot), np.max(x_plot))

    line_centers = []
    line_centers_true = []
    for k in range(len(rb_transitions)):

        peak_ind = np.where( np.abs(x_plot - (rb_transitions[k][0] - cnt_freq + wavemeter_offset)/1e9) < 0.005 )[0][0]
        
        (xf, yf, peak_result, res) = fit_lorentz(x_plot, diff, peak_ind, dn = 5)

        plt.plot(xf, yf, 'k')

        line_centers.append(peak_result.params['cnt0'].value*1e9 + cnt_freq)
        
        # has to have some number of peaks and was given
        line_centers_true.append(rb_transitions[k][0])

        plt.axvline((rb_transitions[k][0] - cnt_freq)/1e9, ls = '--', color = 'r')
        plt.axvline((rb_transitions[k][0] - cnt_freq + wavemeter_offset)/1e9, ls = '--', color = 'k')

    line_centers = np.sort(np.array(line_centers))
    line_centers_true = np.sort(np.array(line_centers_true))

    # plot true value vs fitted line centers)

    plt.subplot(2,1,2) 
    plt.plot((line_centers_true - cnt_freq)/1e9, (line_centers - (line_centers_true + wavemeter_offset))/1e6, 'o-')

    plt.xlabel('True Rb transition frequency (GHz)');
    plt.ylabel('True - Fitted frequency (MHz)');

    plt.xlim(np.min(x_plot), np.max(x_plot))

plt.show()




