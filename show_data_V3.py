""" Original version which saves the csv files
Uses the calculation of finding the frequency spacing between the lines
New modifications: fits the data, and finds params such as temperature, wavemeter offset. Numbers might be off and still needs to be double checked """
"""This version is the current working code, it fits two lines of potassium, the updated fit to include the multiple lines are in the code labeled John, but the fitting is quite off currently """


""" Plan: Plot the fit with the natural line width of the lines"""

""" Current problems:
    color plot no longer works
    
    Questions: 
       1. How does this line work?
       mod_x and mod_y is not passed in, this is the first time they are mentioned
       (mod_x,mod_y) = fcn2min(result.params, fit_x, [], return_plot = True)
       2. What unit is the time in? not sure how that works """
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
main_path = '/Users/kaylajane/software/Prehistoric-Data-Acquisition/'

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

# each slice is delta time units wide, avergaing over the delta time units for each point
def get_time_slot(arr, minx, delta = 20):

    return -np.mean(arr[:, minx:(minx+delta)], axis = 1)

#initialization
res_t = []
res_T = []
res_shift = []


# Time to range over
start_t = 10
end_t   = 100

for k in range(start_t,end_t):

    print(k)
    fit_x = ds * 1e6 + line_act * 1e12
    
    my_t = end_t + k * start_t
    fit_y = get_time_slot(d1, my_t)

    result = my_fit(fit_x, fit_y)
    
    # uses mod_x and mod_y parameters for function to minimize (comes from fit_spectrum)
    (mod_x, mod_y) = fcn2min(result.params, fit_x, [], return_plot = True)

    res_t.append(my_t)
    res_T.append(result.params['T'].value)
    res_shift.append(result.params['shift'].value)

    print('Wavemeter shift = ' + str(result.params['shift'].value/1e6) + 'MHz')
    print('Temperature = ' + str(result.params['T'].value)  + 'K')

    if k == 60:
        plt.figure()
        plt.plot( (fit_x - line_act*1e12)/1e6, fit_y, 'ko',label='data')
        plt.plot( (mod_x - line_act*1e12)/1e6, mod_y, 'r-',label='fitting')
        plt.xlabel('Frequency MHz',fontsize=16)
        plt.ylabel('Strength of line (abs)',fontsize=16)
        plt.tick_params(labelsize=16) #tick size
        plt.tick_params(direction='in') #tick direction
        plt.legend(fontsize=14)
plt.figure()
plt.plot( res_t, res_T, 'mo-') # o is the circle marker
plt.xlabel('Time',fontsize=16)
plt.ylabel('Temperature (K)',fontsize=16)
plt.tick_params(labelsize=16) #tick size
plt.tick_params(direction='in') #tick direction

plt.figure()
plt.plot( res_t, res_shift, 'mo-')
plt.xlabel('Time',fontsize=16)
plt.ylabel('Temperature (K)',fontsize=16)
plt.tick_params(labelsize=16) #tick size
plt.tick_params(direction='in') #tick direction
plt.show()

#plt.xlabel('Frequency Difference (MHz) from from {} GHs'.format(line_act))

plt.figure()

#plt.imshow(d1, aspect = 'auto')
plt.pcolor(dt, ds, d1)

plt.xlabel('Time (ms)')
plt.ylabel('Frequency (MHz)')


plt.show()



d1f.close()
d2f.close()
dsf.close()
