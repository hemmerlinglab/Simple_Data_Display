import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
#from fit_yb import *
#from fit_mo import *
from mpl_toolkits.mplot3d import Axes3D

c = 299792458
yb_174_freq = 751.52653349 # in THz


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

#datafolder = '/Users/boerge/skynet/molecule_computer/'
datafolder = '/home/molecules/software/data/'

#basefolder = str(my_today.strftime('%Y%m%d')) # 20190618
basefolder = '20200518'

basefilename = datafolder + basefolder + '/' + basefolder + '_'

if len(sys.argv)>1:
    time_stamp = sys.argv[1]
else:
    # get latest time stamp
    all_files = np.sort(glob.glob(basefilename + "*"))
    #print(all_files)
    time_stamp = all_files[-1].split('_')[1]


## molybdenum data
#time_stamp = ''


f_posx = basefilename + time_stamp + '_posx'
f_posy = basefilename + time_stamp + '_posy'
f_ch1 = basefilename + time_stamp + '_ch0_arr'


posx = np.genfromtxt(f_posx, delimiter=",")
posy = np.genfromtxt(f_posy, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")


# get number of averages
no_of_avg = 1 # int(len(posx)/len(np.unique(posx)))

print('Found ' + str(no_of_avg) + ' averages.')

posx = av(posx, no_of_avg)
posy = av(posy, no_of_avg)
ch1 = av(ch1, no_of_avg)


inter_x = np.unique(posx)
inter_y = np.unique(posy)



delay_in_for_loop = 100e-6
no_of_time_points = ch1.shape[1]
times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3


# subtracting the DC offset
offset_avg_points = 5
for k in range(ch1.shape[0]):
        #ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -offset_avg_points:-1])
        ch1[k, :] = ch1[k, :] - np.mean(ch1[k, 0:offset_avg_points])

    

cut_time1 = 10.5
cut_time2 = 11.5

ch1_start = np.where( np.abs(times - cut_time1) < 0.5 )[0][0]
ch1_end = np.where( np.abs(times - cut_time2) < 0.5 )[0][0]


target_img = np.zeros([len(inter_x), len(inter_y)])

for nx in range(len(inter_x)):
    for ny in range(len(inter_y)):

        lin_ind = nx * len(inter_y) + ny
        target_img[nx, ny] = np.abs(np.mean(ch1[lin_ind, ch1_start:ch1_end]))


fig = plt.figure(figsize=(10,6))

#plt.subplot(2,1,1)
plt.pcolor(inter_x, inter_y, np.transpose(target_img))
plt.colorbar()

plt.xlabel('x pos')
plt.ylabel('y pos')

#s_inter_x, s_inter_y = np.meshgrid(inter_x, inter_y) 
#
#plt.tight_layout()
#
#plt.subplot(2,1,2)
#
#ax = plt.axes(projection='3d')
#
#ax.plot_surface(s_inter_y, s_inter_x, target_img)
#
#plt.xlabel('x pos')
#plt.ylabel('y pos')








plt.show()

