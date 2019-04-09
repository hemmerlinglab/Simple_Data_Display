import numpy as np
import matplotlib

matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
import csv
import datetime



main_path = '/home/lab-42/data_folder/K_Tests/'

my_today = datetime.datetime.today().strftime('%Y-%m-%d')

my_time = '14-53-13'

f1 = main_path + my_today + '-' + my_time + '_1.csv'
f2 = main_path + my_today + '-' + my_time + '_2.csv'


d1f = open(f1,'r')
d2f = open(f2,'r')

d1r = csv.reader(d1f)
d2r = csv.reader(d2f)

d1 = np.array([])
d2 = np.array([])
for row in d1r:

    hlp_row = ",".join(row).split(',')[:-1]
    hlp = np.array(hlp_row, dtype = np.float)

    if len(hlp)>0:
        d1 = np.vstack((d1, hlp)) if d1.size else hlp

for row in d2r:

    hlp_row = ",".join(row).split(',')[:-1]
    hlp = np.array(hlp_row, dtype = np.float)

    if len(hlp)>0:
        d2 = np.vstack((d2, hlp)) if d2.size else hlp


#d1 = np.array(d1, dtype = np.float)

#plt.figure()
#plt.plot(d1)
#plt.plot(d2)


plt.figure()

a1 = np.mean(d1[10:20, :], axis = 1)
a2 = np.mean(d2[10:20, :], axis = 1)

print(a1)

plt.plot(a1)
plt.plot(a2)

plt.show()



d1f.close()
d2f.close()
