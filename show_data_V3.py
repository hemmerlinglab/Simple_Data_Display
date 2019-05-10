""" Modification to the plots to look much nicer for the presentation """

""" Plan: 
*Abandon transition probability amplitude, too many things we would need to account for

Questions:
	1. does each iteration of k mean 1 ms? Nope

Problems:
	1. cant fix the axes limits without getting indentation errors


"""


import numpy as np
import matplotlib

matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
import csv
import datetime
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import pylab as plb

from fit_spectrum import *

#main_path = '/home/molecules/skynet/Data/K_Tests/'
#main_path = '/home/lab-42/data_folder/K_Tests/'
#main_path = '/home/lab-42/software/github/Prehistoric-Data-Acquisition/'
#main_path = '/home/molecules/software/Prehistoric-Data-Acquisition/'
#main_path = '/Users/johnr/Documents/Github/Prehistoric-Data-Acquisition/'
#main_path = '/Users/kaylajane/software/Prehistoric-Data-Acquisition/'
main_path = 'C:\\Users\\John\\Documents\\Github\\Prehistoric-Data-Acquisition\\'

my_today = datetime.datetime.today().strftime('%Y-%m-%d')

#my_time = '16-32-17'
#my_time = '16-50-54'
#my_time = '17-20-45'

# file names
#f1 = main_path + my_today + '-' + my_time + '_1.csv'
#f2 = main_path + my_today + '-' + my_time + '_2.csv'
#fs = main_path + my_today + '-' + my_time + '_s.csv'
f1 = main_path + 'ch1.txt'
f2 = main_path + 'ch2.txt'
fs = main_path + 'setpts.txt'

# open and read the files
d1f = open(f1,'r')
d2f = open(f2,'r')
dsf = open(fs,'r')

d1r = csv.reader(d1f)
d2r = csv.reader(d2f)
dsr = csv.reader(dsf)

d1 = np.array([])
d2 = np.array([])
ds = np.array([])

for row in d1r:
    # puts data into an array from csv
    hlp_row = ",".join(row).split(',')[:-1]
    hlp = np.array(hlp_row, dtype = np.float)

    # offset correction
    hlp = hlp - hlp[0]


    # vstack -- stack arrays in sequence vertically (row wise)
    if len(hlp)>0:
        d1 = np.vstack((d1, hlp)) if d1.size else hlp

for row in dsr:

	setpoint = ','.join(row).split(',')[0]
	if setpoint != '':
		ds = np.append(ds,float(setpoint))

# ds is the frequency in THz of the wavemeter at the time of any given data point

line_act = 391.01617 #actual frequency of K39 D2 line
ds = ds - line_act # subtract by the actual line

ds = ds * 1e12/1e6 # conversion to MHz

min_volt_ascii = d1.min()

# each slice is delta time units wide, avergaing over the delta time units for each point
def get_time_slot(arr, minx, delta = 20):

    return np.mean(arr[:, minx:(minx+delta)], axis = 1)

def conv_t(ti):
    return np.array([t6*25.0/1500 for t6 in ti])

def conv_A(ab):
    return np.array([(a6+127)*.2*8/256 for a6 in ab])

offset = conv_A([min_volt_ascii])


#
#d1 = np.array(d1, dtype = np.float)
#
#plt.figure()
#plt.plot(d1)
#plt.plot(d2)
#ch2_scale = 2
#d2 = d2*ch2_scale
##
#ch1_scale = .05
#d1 = d1*ch1_scale
dt = np.linspace(0, 2.5, d1.shape[1]) * 10.0
#
# Oscilloscope Figure
plt.figure()
#
plt.plot(dt, conv_A(d1[25, :]))

plt.xlabel('Time (ms)',fontsize=16)
plt.ylabel('Voltage',fontsize=16)
plt.tick_params(labelsize=16) #tick size
plt.tick_params(direction='in') #tick direction
plt.tight_layout()
## figure 2 spectrum stuff
plt.figure()
#

#initialization
res_t = []
res_T = []
res_shift = []



# Time to range over
start_t = 10
end_t   = 100
tick_array = np.linspace(-600,600,num=7)

for k in range(start_t,end_t):

    print(k)
    fit_x = ds * 1e6 + line_act * 1e12
    
    # time in ms
    my_t = end_t + k * start_t
    fit_y = get_time_slot(d1, my_t)

    result = my_fit(fit_x, fit_y)
    
    # uses mod_x and mod_y parameters for function to minimize (comes from fit_spectrum)
    # defines variables mod_x and mod_y, and fills with the return of that function
    (mod_x, mod_y) = fcn2min(result.params, fit_x, [], return_plot = True)

    res_t.append(my_t)
    res_T.append(result.params['T'].value)
    res_shift.append(result.params['shift'].value)

    print('Wavemeter shift = ' + str(result.params['shift'].value/1e6) + 'MHz')
    print('Temperature = ' + str(result.params['T'].value)  + 'K')


    if k == 60:
        plt.figure()

        fit_y_1 = conv_A(fit_y) - offset[0]
        mod_y_1 = conv_A(mod_y) - offset[0]
        plt.plot( (fit_x - line_act*1e12)/1e6, fit_y_1, 'ko',label='data')
        plt.plot( (mod_x - line_act*1e12)/1e6, mod_y_1, 'r-',label='fitting')
        plt.xlabel('Frequency MHz',fontsize=16)
        plt.ylabel('Precent of Signal Absorbed (abs)',fontsize=16)
        plt.tick_params(labelsize=16) #tick size
        plt.tick_params(direction='in') #tick direction
        plt.xlim(-700,700)
        plt.xticks(tick_array)
        plt.legend(fontsize=14)

print(len(res_t))
print(len(res_T))

# only take reliable data
# after a certain amount of time, we are off target so does not make sense to plot that temperature
strt_cut = 12
end_cut  = 50
res_t = conv_t(res_t)

plt.figure()
plt.plot(res_t[strt_cut:end_cut], res_T[strt_cut:end_cut], 'mo-') # o is the circle marker
plt.xlabel('Time (ms)',fontsize=16)
plt.ylabel('Temperature (K)',fontsize=16)
plt.tick_params(labelsize=16) #tick size
plt.tick_params(direction='in') #tick direction

plt.figure()
plt.plot( res_t, res_shift, 'mo-')
plt.xlabel('Time (ms)',fontsize=16)
plt.ylabel('Temperature (K)',fontsize=16)
plt.tick_params(labelsize=16) #tick size
plt.tick_params(direction='in') #tick direction
plt.show()

#plt.xlabel('Frequency Difference (MHz) from from {} GHs'.format(line_act))

plt.figure()

#plt.imshow(d1, aspect = 'auto')
plt.pcolor(dt, ds, d1)

plt.xlabel('Time (ms)',fontsize=16)
plt.ylabel('Frequency (MHz)',fontsize=16)
plt.tick_params(labelsize=16) #tick size
plt.tick_params(direction='in') #tick direction
plt.yticks(tick_array)

plt.tight_layout()

plt.show()



d1f.close()
d2f.close()
dsf.close()
