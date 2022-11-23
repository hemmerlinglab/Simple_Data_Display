import numpy as np
import matplotlib.pyplot as plt

import sys

basefilename = '/home/electrons/software/data/20221123/20221123_'

filename = basefilename + sys.argv[1] + '_scan_result'

d = np.genfromtxt(filename, delimiter = ',')

# change units to seconds
d = d*1e-6

filename = basefilename + sys.argv[1] + '_set_points'

f = np.genfromtxt(filename, delimiter = ',')

print(f)
f = (f - np.mean(f))*1e12/1e6

no_of_bins = 30
det_time = 20e-3

t = np.linspace(0, det_time, no_of_bins)


my_histograms = np.zeros([len(f), len(t)])
for k in range(d.shape[0]):
    (hlp, my_bins) = np.histogram(d[k, :], bins = t)

    my_histograms[k, :-1] = hlp

# remove the zero bin
my_histograms = my_histograms[:, 1:]

#plt.pcolor(t[1:]/1e-3, f, my_histograms[1:, :])
plt.pcolor(t/1e-3, f, my_histograms)


plt.xlabel('Time (ms)')

plt.ylabel('Frequency (MHz)')

plt.show()


