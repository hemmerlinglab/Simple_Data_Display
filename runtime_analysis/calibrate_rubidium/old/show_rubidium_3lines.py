import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys

from fit_rb_3lines import *
from scipy.signal import find_peaks

c = 299792458


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

#datafolder = '/Users/boerge/software/data/molecule_computer/'
#datafolder = '/home/molecules/software/data/'
datafolder = '/Users/boerge/Software/offline_data/'

basefolder = '20200303'

time_stamp = sys.argv[1] #'145643'

basefilename = datafolder + basefolder + '/' + basefolder + '_'


f_freqs = basefilename + time_stamp + '_set_points'
freqs = np.genfromtxt(f_freqs, delimiter=",")
# get number of averages
no_of_avg = int(len(freqs)/len(np.unique(freqs)))
print('Found ' + str(no_of_avg) + ' averages.')

for k in range(3):
    exec('f_ch' + str(k) + ' = "' + basefilename + time_stamp + '_ch' + str(k) + '_arr"')
    exec('ch' + str(k) + ' = np.genfromtxt(f_ch' + str(k) + ', delimiter=",")')
    exec('ch' + str(k) + ' = av(ch' + str(k) + ', no_of_avg)')

freqs = av(freqs, no_of_avg)

ch0_mean = np.mean(ch0, axis = 1)
ch1_mean = np.mean(ch1, axis = 1)
ch2_mean = np.mean(ch2, axis = 1)

ch0_mean -= ch0_mean[0]
ch1_mean -= ch1_mean[0]
ch2_mean -= ch2_mean[0]


# switch 0 and 1
#ch0_mean = ch1_mean



diff_sig = ch0_mean - ch1_mean


delay_in_for_loop = 100e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3

   
Rb87 = {
        'atom' : 'Rb87',
        'S1/2-P3/2' : 384.230484468e12,
        'S1/2' : {
            '2' : -2.563005e9,
            '1' : 4.27167e9
            },
        'P3/2' : {
            '3' : +193.7407e6,
            '2' : -72.9112e6,
            '1' : -229.8518e6,
            '0' : -302.0738e6
            }
    }

Rb85 = {
        'atom' : 'Rb85',
        'S1/2-P3/2' : 384.230406373e12,
        'S1/2' : {
            '3' : -1.264888516e9,
            '2' : 1.770843922e9
            },
        'P3/2' : {
            '4' : +100.205e6,
            '3' : -20.435e6,
            '2' : -83.835e6,
            '1' : -113.208e6
            }
    }


def create_transition(atom, gs, es, Fg, Fe, Fe_co = '', crossover = False):

    if not crossover:
        freq = atom[gs + '-' + es] + atom[gs][str(Fg)] + atom[es][str(Fe)]
        label = atom['atom'] + ' ' + Fg + '->' + Fe
        ls = '-'
    else:
        freq = atom[gs + '-' + es] + atom[gs][str(Fg)] + (atom[es][str(Fe)] + atom[es][str(Fe_co)])/2.0
        label = atom['atom'] + ' ' + Fg + '->' + Fe + '/' + Fe_co + ' (co)'
        ls = '--'

    return [freq, label, ls]


setpoint_offset = 384.230 # should be read out from the config file


# only keeping the lines that cycle photons

my_lines = []
my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '3'))
#my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '2'))
#my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '1'))

#my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '1', Fe_co = '2', crossover = True))
my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '2', Fe_co = '3', crossover = True))
my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '1', Fe_co = '3', crossover = True))

my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '4'))
#my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '3'))
#my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '2'))

#my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '2', Fe_co = '3', crossover = True))
my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '3', Fe_co = '4', crossover = True))
my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '2', Fe_co = '4', crossover = True))


my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '3'))
##my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '2'))
##my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '1'))

#my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '1', Fe_co = '2', crossover = True))

my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '2', Fe_co = '3', crossover = True))
my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '1', Fe_co = '3', crossover = True))

# fit spectrum

(xfit, yfit, fit_result) = fit_rb(freqs, ch0_mean, my_lines, setpoint_offset)
#(xfit, yfit, fit_result) = fit_rb(freqs[0:200], ch0_mean[0:200], my_lines, setpoint_offset)

wavemeter_offset = fit_result.params['x_offset'].value
print('Wavemeter offset: {}'.format(wavemeter_offset))




myfontsize = 16


ind = []
for k in range(len(my_lines)):
    ind.append( (my_lines[k][0]/1.0e12 - setpoint_offset )*1.0e6 - wavemeter_offset )

mymax = np.max(diff_sig)

plt.figure(figsize=(10,6))
#plt.subplot(2,1,1)
#plt.plot(freqs, ch1_mean)
#plt.plot(freqs, ch2_mean)

plt.xlim(np.min(freqs), np.max(freqs))

plt.xlabel('Freqs (MHz) + ' + str(setpoint_offset) + ' THz', fontsize = myfontsize)

plt.ylabel('Absorption Signal (a.u)', fontsize = myfontsize)
plt.tick_params(labelsize=14, direction='in')


dv = 1.0/len(my_lines)
for k in range(len(my_lines)):
    plt.axvline(ind[k], ls = my_lines[k][2], color = 'r')
    plt.text(ind[k] + 10, 0.0, my_lines[k][1], rotation = 90, fontsize = 8)

plt.plot(freqs, ch0_mean, linewidth = 2)

#plt.subplot(2,1,2)
#plt.figure()
#plt.plot(freqs, diff_sig)
#
#plt.xlabel('Freqs (MHz) + ' + str(setpoint_offset) + ' THz')
#
#plt.xlim(np.min(freqs), np.max(freqs))


plt.plot(xfit, yfit, 'r')




plt.show()




