""" Original version which saves the csv files
Uses the calculation of finding the frequency spacing between the lines
New modifications: fits the data, and finds params such as temperature, wavemeter offset. Numbers might be off and still needs to be double checked """
"""This version is the current working code, it fits two lines of potassium, the updated fit to include the multiple lines are in the code labeled John, but the fitting is quite off currently """


""" Plan: Plot the fit with the natural line width of the lines"""
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
main_path = '/home/molecules/software/Prehistoric-Data-Acquisition/'
#main_path = '/Users/johnr/Documents/Github/Prehistoric-Data-Acquisition/'

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

#for row in d2r:
#
#    hlp_row = ",".join(row).split(',')[:-1]
#    hlp = np.array(hlp_row, dtype = np.float)
#
#    if len(hlp)>0:
#        print(len(hlp))
#        d2 = np.vstack((d2, hlp)) if d2.size else hlp

for row in dsr:

	setpoint = ','.join(row).split(',')[0]
	if setpoint != '':
		ds = np.append(ds,float(setpoint))

line_act = 391.01617
ds = ds - line_act

ds = ds * 1e12/1e6

def get_time_slot(arr, minx, delta = 20):

    return -np.mean(arr[:, minx:(minx+delta)], axis = 1)


res_t = []
res_T = []
res_shift = []

# fitting
for k in range(10,100):

    print(k)
    fit_x = ds * 1e6 + line_act * 1e12
    
    my_t = 100 + k * 10
    fit_y = get_time_slot(d1, my_t)

    result = my_fit(fit_x, fit_y)

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
        plt.ylabel('Strength of line (arb)',fontsize=16)
        plt.tick_params(labelsize=16) #tick size
        plt.tick_params(direction='in') #tick direction
        plt.legend(fontsize=14)
plt.figure()
plt.plot( res_t, res_T, 'mo-') # o is the circle marker

plt.xlabel('Time',fontsize=16)
plt.ylabel('Temperature (K)',fontsize=16)
plt.tick_params(labelsize=16) #tick size
plt.tick_params(direction='in') #tick direction
plt.legend(fontsize=14)
plt.figure()
plt.plot( res_t, res_shift, 'o-')

plt.show()


