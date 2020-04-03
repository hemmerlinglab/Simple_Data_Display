import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
from helper import *
from fit_rb import *


my_date = '20200303'
my_time = '121749'

(freqs, signal) = get_data(my_date, my_time, datafolder = '/Users/boerge/Software/offline_data/')


rb_transitions = calculate_Rb_transitions()

# fit spectrum

(xfit, yfit, fit_result) = fit_rb(freqs, signal[1], rb_transitions)

wavemeter_offset = fit_result.params['x_offset'].value
print('Wavemeter offset: {0:2.6} MHz'.format(wavemeter_offset/1e6))




cnt_freq = 384.230484468e12 # S1/2 - P3/2
myfontsize = 16

plt.figure(figsize = (20,6))
plt.plot((freqs - cnt_freq)/1e9, signal[1], 'r-')
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



plt.title('Rb Calibration Scan - ' + my_date + '_' + my_time)

plt.tight_layout()

plt.savefig(my_date + '_' + my_time + '_Rb_calibration.png', format = 'png')

plt.show()




