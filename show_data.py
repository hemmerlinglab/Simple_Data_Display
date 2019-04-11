import numpy as np
import matplotlib

matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
import csv
import datetime


main_path = '/home/molecules/skynet/Data/K_Tests/'
#main_path = '/home/lab-42/data_folder/K_Tests/'

my_today = datetime.datetime.today().strftime('%Y-%m-%d')

my_time = '17-27-52'

f1 = main_path + my_today + '-' + my_time + '_1.csv'
f2 = main_path + my_today + '-' + my_time + '_2.csv'
fs = main_path + my_today + '-' + my_time + '_set.csv'


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

    hlp_row = ",".join(row).split(',')[:-1]
    hlp = np.array(hlp_row, dtype = np.float)

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


ds = ds - 391.0160

ds = ds * 1e12/1e6

#d1 = np.array(d1, dtype = np.float)

#plt.figure()
#plt.plot(d1)
#plt.plot(d2)
#ch2_scale = 2
#d2 = d2*ch2_scale

ch1_scale = .05
d1 = d1*ch1_scale



dt = np.linspace(0, 2.5, d1.shape[1]) * 10.0

print(dt)

plt.figure()

plt.plot(dt, d1[4, :])

plt.figure()

minx = 510
maxx = 520

a1 = np.mean(d1[:, minx:maxx], axis = 1)

plt.plot(ds, a1)

plt.xlabel('Frequency (MHz)')

plt.figure()

#plt.imshow(d1, aspect = 'auto')
plt.pcolor(dt, ds, d1)

plt.xlabel('Time (ms)')
plt.ylabel('Frequency (MHz)')


plt.show()



d1f.close()
d2f.close()
dsf.close()
