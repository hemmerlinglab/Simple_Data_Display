import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
#from fit_yb import *
#from fit_mo import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import *
from helper import *



def integrate_absorption(img, inter_x, inter_y, times, time_cut):

    dt = times[1] - times[0]
    t_start = np.where( np.abs(times - time_cut[0]) < dt )[0][0]
    t_stop  = np.where( np.abs(times - time_cut[1]) < dt )[0][0]
   
    abs_img = np.zeros([len(inter_x), len(inter_y)])

    for nx in range(len(inter_x)):
        for ny in range(len(inter_y)):
    
            lin_ind = nx * len(inter_y) + ny
            
            absorption = np.mean(img[lin_ind, t_start:t_stop])
            
            abs_img[nx, ny] = np.abs(absorption)

    return np.transpose(abs_img)



def plot_single_image(inter_x, inter_y, img, color_max = 1.0, factor = 1000.0, title = ''):

    plt.pcolor(inter_x, inter_y, factor * img)
    plt.colorbar()
    plt.clim(0, factor * color_max)
    
    plt.xlabel('x pos')
    plt.ylabel('y pos')

    plt.gca().invert_yaxis()
    
    plt.title(title)

    #plt.tight_layout()

    return

#################################################################
# main
#################################################################

(inter_x, inter_y, times, ch1) = get_data(sys.argv[1], sys.argv[2])


#basefolder = '20200724'
#time_stamp = '155526'

# get images for foreground and background

t1 = 10.0
t2 = 11.0

t1_bg = 30.0
t2_bg = 39.0

target_img = integrate_absorption(ch1, inter_x, inter_y, times, [t1, t2])

bg_img = integrate_absorption(ch1, inter_x, inter_y, times, [t1_bg, t2_bg])



color_max = np.max(np.max(target_img))

filtered_img = uniform_filter(target_img, size=2, mode='constant')





fig = plt.figure(figsize=(10,10))

plt.subplot(2,1,1)

plot_single_image(inter_x, inter_y, target_img, color_max = color_max, title = "Target, t = {0:.1f} - {1:.1f} ms".format(t1, t2))

plt.subplot(2,1,2)

plot_single_image(inter_x, inter_y, bg_img, color_max = color_max, title = "Background, t = {0:.1f} - {1:.1f} ms".format(t1_bg, t2_bg))



plt.show()







