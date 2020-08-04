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


def get_avg_absorption(coords,super_image):
	# coords is in index format with coords = [ xmin, xmax, ymin, ymax ]
	xmin = coords[2]
	xmax = coords[3]
	ymin = coords[0]
	ymax = coords[1]
	sub_image = target_img[xmin:xmax,ymin:ymax]
	# print(sub_image.shape)
	return np.mean(sub_image),sub_image


#################################################################
# main
#################################################################

(inter_x, inter_y, times, ch1, conf) = get_img_data(sys.argv[1], sys.argv[2])

#basefolder = '20200724'
#time_stamp = '155526'

# get images for foreground and background

t1 = 10.0
t2 = 11.0

t1_bg = 30.0
t2_bg = 39.0

target_img = integrate_absorption(ch1, inter_x, inter_y, times, [t1, t2])

bg_img = integrate_absorption(ch1, inter_x, inter_y, times, [t1_bg, t2_bg])

color_max = np.max(np.max(target_img)) * 0.7

# print(target_img.shape)
c1_vals = [2, 10, 10, 24]
mean1,sub_im1 = get_avg_absorption(c1_vals,target_img)
c2_vals = [12, 22, 8, 28]
mean2,sub_im2 = get_avg_absorption(c2_vals,target_img)

# Left: 0.05g Al, 0.45g KCl
# Right: 0.25g AlCl3

# Molar weights:
Al = 26.982
K = 39.098
Cl = 35.45

Al_KCl = np.min((0.05/Al, 0.45/(K + Cl)))/0.5
AlCl3 = 1/(Al + 3 * Cl)

print('Al+KCl: ',Al_KCl,'mols, AlCl3: ',AlCl3,'mols')
print('Means: ',mean1,mean2)


abs_comp1 = mean1/Al_KCl
abs_comp2 = mean2/AlCl3
print('Absorption Factors: ',abs_comp1,abs_comp2)
abs_comp = np.max((abs_comp1/abs_comp2,abs_comp2/abs_comp1))

if abs_comp1 > abs_comp2:
	print('Left is {} times greater than Right'.format(str(abs_comp)))
elif abs_comp2 > abs_comp1:
	print('Right is {} times greater than Left'.format(str(abs_comp)))
else:
	print('Identical absorption... illogical...')



# print(inter_y.shape)
# print('Left: {}, Right: {}'.format(str(mean1),str(mean2)))
# plt.figure()
# plt.subplot(1,2,1)
# plot_single_image(inter_x[c1_vals[0]:c1_vals[1]],inter_y[c1_vals[2]:c1_vals[3]],sub_im1,color_max=color_max,title='Left Target')
# plt.subplot(1,2,2)
# plot_single_image(inter_x[c2_vals[0]:c2_vals[1]],inter_y[c2_vals[2]:c2_vals[3]],sub_im2,color_max=color_max,title='Right Target')

filtered_img = uniform_filter(target_img, size=2, mode='constant')



fig = plt.figure(figsize=(10,3))

plt.subplot(1,2,1)

plot_single_image(inter_x, inter_y, target_img, color_max = color_max, title = "Target, t = {0:.1f} - {1:.1f} ms".format(t1, t2))

x0 = 5.15
y0 = 4.95
r = 0.55

tpar = np.pi/180 * np.linspace(40, 320, 100)

x = r * np.cos(tpar) + x0
y = r * np.sin(tpar) + y0

plt.plot(x, y, 'r--')
plt.subplot(1,2,2)

plot_single_image(inter_x, inter_y, bg_img, color_max = color_max, title = "Background, t = {0:.1f} - {1:.1f} ms".format(t1_bg, t2_bg))




fig = plt.figure(figsize=(15,3))

t_arr = np.linspace(0,4.0,5)

for dt,n in enumerate(t_arr):
    plt.subplot(1,5,n+1)

    t1 = 10.0
    t2 = 10.1 + dt
    target_img = integrate_absorption(ch1, inter_x, inter_y, times, [t1, t2])

    plot_single_image(inter_x, inter_y, target_img, color_max = color_max, title = "Target, t = {0:.1f} - {1:.1f} ms".format(t1, t2))

plt.show()







