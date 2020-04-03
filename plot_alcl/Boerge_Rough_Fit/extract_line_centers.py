import numpy as np
import pickle
import matplotlib.pyplot as plt
from fit_func import *
from helper import *




# load data

arr = pickle.load( open( "data.pickle", "rb" ) )

f_arr = arr['f_arr']
sig_arr = arr['sig_arr']
data = arr['data']


c = 299792458

cnt_freq = 100.0 * c * 38237.00


# fit each data set with a gaussian to find center of lines (roughly)

line_pos = []

x_data = []
offset_free_data = []

######################################
# P lines (v -> v') = 0->0
######################################

plt.figure()
for n in [0,1,2,3,4,5]:

    x = f_arr[n]
    y = sig_arr[n]

    x_mean = np.mean(x)
    x_min = np.min(x)
    x_max = np.max(x)
    y_min = np.min(y)
    y_max = np.max(y)

    params = Parameters()

    params.add('p0', value=y_min, min=-1.0, max=1.0, vary = True)
    params.add('p1', value=(y_max - y_min), min=0.0, max=3.0, vary = True)
    params.add('p2', value=x_mean, min=x_min-10e9, max=x_mean+10e9, vary = True)
    params.add('p3', value=(x_max-x_min)/10.0, min=0.001, max=(x_max-x_min), vary = True)
 
    (xf, yf, m) = fit_func(x, y, params, gauss)

    line_pos.append(m.params['p2'].value)

    x_data.append(x)
    offset_free_data.append((y - m.params['p0'])/m.params['p1'])

    plt.subplot(3,3,n+1)
    plt.plot((x - cnt_freq)/1e9, y, 'ro', label = data[n]['time'])
    plt.plot((xf - cnt_freq)/1e9, yf, 'k-')

    plt.xlabel("f (GHz) - {0:2.3f} THz".format(cnt_freq/1e12))
    plt.legend()
    plt.tight_layout()


#######################################################################
# Q line(s) (v -> v') = 0->0, need to be fitted to a sum of gaussians
#######################################################################


x = []
y = []
for n in [6,7]:

    x.extend(f_arr[n])
    y.extend(sig_arr[n])

x = np.array(x)
y = np.array(y)

x_mean = np.mean(x)
x_min = np.min(x)
x_max = np.max(x)
y_min = np.min(y)
y_max = np.max(y)

params = Parameters()

params.add('p0', value=y_min, min=-1.0, max=1.0, vary = True)
params.add('p3', value=(x_max-x_min)/10.0, min=0.001, max=(x_max-x_min), vary = True)

for k in range(4):
    # amplitudes
    params.add('p1' + str(k), value=(y_max - y_min), min=0.0, max=3.0, vary = True)

    # positions
    params.add('p2' + str(k), value=x_min + 1.0e9 * k, min=x_min-10e9, max=x_mean+10e9, vary = True)


(xf, yf, m) = fit_func(x, y, params, multi_gauss)

for k in range(4):
    line_pos.append(m.params['p2' + str(k)].value)

x_data.append(x)
offset_free_data.append( (y - m.params['p0']) / (m.params['p10']+m.params['p11']) )

plt.figure()
plt.plot((x - cnt_freq)/1e9, y, 'ro')
plt.plot((xf - cnt_freq)/1e9, yf, 'k-')

plt.xlabel("f (GHz) - {0:2.3f} THz".format(cnt_freq/1e12))
plt.tight_layout()




###########################################################################
# Q line(s) (v -> v') = 0->0, need to be fitted to a sum of gaussians
#########################################################3#################

line_pos = np.sort(np.array(line_pos))

#plt.figure()
#plt.plot((line_pos - cnt_freq)/1e9, 'x')





##################################################################################################
# get all fitted lines and subtract offset, then compare to simulated data by varying Be and Te
##################################################################################################

# take only unique values
uni_lines = line_pos[[0,1,2,5,6]]


gauss_data = []

ampl = 1.0
width = 0.5e9
gauss_x = cnt_freq + np.linspace(-100e9, 30e9, 500)

for k in range(len(uni_lines)):

    gauss_data.append(simple_gauss(gauss_x, uni_lines[k], ampl, width))

gauss_data = np.array(gauss_data)

x_data = np.array(x_data)
offset_free_data = np.array(offset_free_data)





with open('line_centers.pickle', 'wb') as f:
    pickle.dump({'line_pos' : line_pos, 'unique_lines' : uni_lines, 'x_data' : x_data, 'offset_free_data' : offset_free_data}, f)



plt.show()


