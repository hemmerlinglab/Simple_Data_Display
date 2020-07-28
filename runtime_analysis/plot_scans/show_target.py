import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
#from fit_yb import *
#from fit_mo import *
from mpl_toolkits.mplot3d import Axes3D
from configparser import *

c = 299792458
yb_174_freq = 751.52653349 # in THz


from scipy.ndimage import *

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


def read_in_config(f):
    
    config = ConfigParser()
    config.read(f)

    sensor_ids = config.sections()
    # make dictionary out of config

    sensors = {}

    for s in sensor_ids:
        opts = config.options(s)
        
        sensors[s] = {}
        for o in opts:
            sensors[s][o] = config.get(s, o)

    return sensors




my_today = datetime.datetime.today()

datafolder = '/Users/boerge/software/offline_data/'
#datafolder = '/home/molecules/software/data/'

#basefolder = str(my_today.strftime('%Y%m%d')) # 20190618
basefolder = '20200728'
time_stamp = '141306'

#basefolder = '20200724'
#time_stamp = '155526'

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
f_conf = basefilename + time_stamp + '_conf'


posx = np.genfromtxt(f_posx, delimiter=",")
posy = np.genfromtxt(f_posy, delimiter=",")
ch1 = np.genfromtxt(f_ch1, delimiter=",")

conf = read_in_config(f_conf)

# get number of averages
no_of_avg = int(conf['scan_count']['val'])
print('Found ' + str(no_of_avg) + ' averages.')


#posx = av(posx, no_of_avg)
#posy = av(posy, no_of_avg)
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

    

cut_time1 = 10.0
cut_time2 = 10.5

shift = 20.0
cut_time1_offset = 10.0+shift
cut_time2_offset = 10.5+shift


ch1_start = np.where( np.abs(times - cut_time1) < 0.25 )[0][0]
ch1_end = np.where( np.abs(times - cut_time2) < 0.25 )[0][0]


ch1_start_offset = np.where( np.abs(times - cut_time1_offset) < 0.25 )[0][0]
ch1_end_offset = np.where( np.abs(times - cut_time2_offset) < 0.25 )[0][0]



target_img = np.zeros([len(inter_x), len(inter_y)])

bg_img = np.zeros([len(inter_x), len(inter_y)])



for nx in range(len(inter_x)):
    for ny in range(len(inter_y)):

        lin_ind = nx * len(inter_y) + ny
        
        absorption = np.mean(ch1[lin_ind, ch1_start:ch1_end])
        
        bg_signal = np.mean(ch1[lin_ind, ch1_start_offset:ch1_end_offset])

        target_img[nx, ny] = np.abs(absorption)
        
        bg_img[nx, ny] = np.abs(bg_signal)

        #plt.plot(times, ch1[lin_ind, :])
        #plt.axvline(times[ch1_start])
        #plt.axvline(times[ch1_end])
        #plt.show()
        #asd





color_max = np.max(np.max(target_img))

filtered_img = uniform_filter(target_img, size=5, mode='constant')


filtered_img = target_img - bg_img



fig = plt.figure(figsize=(10,7))

plt.subplot(2,2,1)
#plt.imshow(np.transpose(target_img), cmap='Blues', interpolation='nearest')
plt.pcolor(inter_x, inter_y, np.transpose(target_img))
plt.colorbar()

plt.clim(0, color_max)

plt.xlabel('x pos')
plt.ylabel('y pos')

plt.title('High value = High absorption')


plt.gca().invert_yaxis()


plt.subplot(2,2,2)
#plt.imshow(np.transpose(bg_img), cmap='Blues', interpolation='nearest')
plt.pcolor(inter_x, inter_y, np.transpose(bg_img))
plt.colorbar()

plt.clim(0, color_max)

plt.xlabel('x pos')
plt.ylabel('y pos')

plt.title('High value = High absorption')


plt.gca().invert_yaxis()



#fig = plt.figure(figsize=(10,6))

plt.subplot(2,2,3)
#plt.imshow(np.transpose(filtered_img), cmap='Blues', interpolation='nearest')
plt.pcolor(inter_x, inter_y, np.transpose(filtered_img))
plt.colorbar()

plt.clim(0, color_max)

plt.xlabel('x pos')
plt.ylabel('y pos')

plt.title('High value = High absorption')

plt.gca().invert_yaxis()

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


#plt.figure()
#plt.pcolor(inter_x, inter_y[::3], np.transpose(target_img[:,::3]))

plt.figure()

hlp = np.transpose(target_img)
y_sum = np.sum(hlp, axis = 1)
x_sum = np.sum(hlp, axis = 0)

plt.subplot(2,1,1)
plt.plot(x_sum)
plt.subplot(2,1,2)
plt.plot(y_sum)

fig = plt.figure(figsize=(10,7))

plt.subplot(2,1,1)
plt.imshow(np.transpose(target_img), cmap='Blues', interpolation='nearest')

plt.subplot(2,1,2)
plt.imshow(np.transpose(filtered_img), cmap='Blues', interpolation='nearest')



plt.figure()

y = np.mean(ch1, axis = 0)



plt.plot(times, y)


plt.show()







