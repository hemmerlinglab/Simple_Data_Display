""" Original version which saves the csv files
Uses the calculation of finding the frequency spacing between the lines
New modifications: fits the data, and finds params such as temperature, wavemeter offset. Numbers might be off and still needs to be double checked """

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
main_path = '/home/lab-42/software/github/Prehistoric-Data-Acquisition/'
#main_path = '/home/molecules/software/Prehistoric-Data-Acquisition/'
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

#d1 = np.array(d1, dtype = np.float)

#plt.figure()
#plt.plot(d1)
##plt.plot(d2)
##ch2_scale = 2
##d2 = d2*ch2_scale
#
#ch1_scale = .05
#d1 = d1*ch1_scale
#dt = np.linspace(0, 2.5, d1.shape[1]) * 10.0
#
#plt.figure()
#
#plt.plot(dt, d1[25, :])
#
## figure 2 spectrum stuff
#plt.figure()
#

def get_time_slot(arr, minx, delta = 20):

    return -np.mean(arr[:, minx:(minx+delta)], axis = 1)


res_t = []
res_T = []
res_shift = []

# fitting


#params = Parameters()
#params.add('a1', value=1.0, min = 0)
#params.add('a2', value=1.0, min = 0)
#params.add('a3', value=1.0, min = 0)
#params.add('a4', value=1.0, min = 0)
#params.add('a5', value=1.0, min = 0)
#params.add('shift', value=-5.0e6, min = -500.0e6, max = 500.0e6)
#params.add('T', value=2.0, min = 0.1, max = 50.0)
#
#fit_x = ds * 1e6 + line_act * 1e12
#
#fit_x = line_act * 1e12 + np.linspace(-200, 600, 100) * 1e6
#(mod_x, mod_y) = fcn2min(params, fit_x, [], return_plot = True)
#
#plt.figure()
#plt.plot((mod_x - line_act*1e12)/1e6, mod_y)
#plt.show()
#asd

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
        plt.plot( (fit_x - line_act*1e12)/1e6, fit_y, 'ko')
        plt.plot( (mod_x - line_act*1e12)/1e6, mod_y, 'r-')
        plt.xlabel('Frequency MHz',fontsize=16)
        plt.ylabel('Strength of line (arb)',fontsize=16)
plt.figure()
plt.plot( res_t, res_T, 'o-')

plt.figure()
plt.plot( res_t, res_shift, 'o-')

plt.show()

asd


# plot the natural line width against plot we have here

peak1x = int((len(ds)/2))
""" Fit to a gaussian """
def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

popt,pcov = curve_fit(gaus,ds[0:peak1x],arb[0:peak1x])

plt.plot(ds[0:peak1x], gaus(ds[0:peak1x],*popt),'r')

non_hyper39 = 766.701
pf3sf2_39 = 14.4-173.1
pf2sf1_39 = 288.6-6.7
freq_diff = abs(pf3sf2_39-pf2sf1_39)
print('Frequency Difference MHz: {}'.format(freq_diff))

#
#lines = [pf3sf2_39,pf2sf1_39]
#for line in lines:
#	plt.axvline(x=line)
#

plt.xlabel('Frequency Difference (MHz) from from {} GHs'.format(line_act))

plt.figure()

#plt.imshow(d1, aspect = 'auto')
plt.pcolor(dt, ds, d1)

plt.xlabel('Time (ms)')
plt.ylabel('Frequency (MHz)')


plt.show()



d1f.close()
d2f.close()
dsf.close()
