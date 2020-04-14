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



wavemeter_offset = -32.6e6

# apply wavemeter offset to all x_values
for k in range(len(f_arr)):

    f_arr[k] = f_arr[k] - wavemeter_offset



c = 299792458

cnt_freq = 100.0 * c * 38237.00# + 15.5e9





# fit each data set with a Gaussian to find center of lines (roughly)

line_pos = []

x_data = []
offset_free_data = []

xf_data = []
yf_data = []

######################################
# P lines (v -> v') = 0->0
######################################

#plt.figure()
#for n in [0,1,2,3,4,5]:
for n in [2,3,4,5]:

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

    xf_data.append(xf)
    yf_data.append( (yf - m.params['p0'])/m.params['p1'] )

    #plt.subplot(3,3,n+1)
    #plt.plot((x - cnt_freq)/1e9, y, 'ro', label = data[n]['time'])
    #plt.plot((xf - cnt_freq)/1e9, yf, 'k-')

    #plt.xlabel("f (GHz) - {0:2.3f} THz".format(cnt_freq/1e12))
    #plt.legend()
    #plt.tight_layout()




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

xf_data.append(xf)
yf_data.append( (yf - m.params['p0']) / (m.params['p10']+m.params['p11']+m.params['p12']+m.params['p13']) * 2.0)

hlp = []
for k in range(4):
    hlp.append(m.params['p2' + str(k)].value)

line_pos.append(min(hlp))

x_data.append(x)
offset_free_data.append( (y - m.params['p0']) / (m.params['p10']+m.params['p11']+m.params['p12']+m.params['p13']) * 2.0)





# set cnt_freq for plotting to Q00 line

cnt_freq = hlp[1]





from brokenaxes import brokenaxes

plt.figure()

plot_x = (np.array(x_data) - cnt_freq)/1e9

plot_x = [l.tolist() for l in plot_x]


broken_xlims = []
for k in range(len(plot_x)):
    
    avg_cnt = (plot_x[k][0] + plot_x[k][-1])/2.0
    offset = 1.5
    hlp = [avg_cnt - offset, avg_cnt + offset]
    broken_xlims.append(hlp)

broken_xlims.sort(key = lambda x : x[0])



line_pos = np.sort(line_pos)

my_color = ['r', 'k', 'b', 'g', 'y']

bax = brokenaxes(xlims=broken_xlims) #, subplot_spec=sps1)

for k in range(len(x_data)):

    
    bax.scatter(plot_x[k], offset_free_data[k], color = my_color[k], marker = 'o', fc = 'w', lw = 2.0)

    bax.plot(plot_x[k], offset_free_data[k], color = my_color[k], marker = 'o', ls = '', alpha = 0.3, mew = 0.0)
    
    bax.plot( (xf_data[k] - cnt_freq)/1e9 , yf_data[k], '-', color = my_color[k], label = line_pos[k]/1e12)

    
    #bax.text( (line_pos[0] - cnt_freq)/1e9, 1.2, 'asd' )
    
    bax.axvline( (line_pos[k] - cnt_freq)/1e9 , ls = '--', ymax = 0.1)
    
    #bax.text( (line_pos[k] - cnt_freq)/1e9, -0.1, "{0:2.6f}".format(line_pos[k]/1e12) )


plt.text( -29, 1.0, "{0:2.6f}".format(line_pos[k]/1e12) )

bax.set_xlabel("Frequency (GHz) + {0:2.6f} THz".format(cnt_freq/1e12), labelpad = 20)
bax.set_ylabel("Absorption Signal (a.u., normalized)", labelpad = 25)

for k in range(len(line_pos)):
    print("{0:2.6f} THz".format(line_pos[k]/1e12))



#plt.figure()
#
#for n in range(len(f_arr)):
#    plt.plot(f_arr[n], sig_arr[n])


plt.show()


# save data in ascii files



for k in range(len(line_pos)):
    np.savetxt('x' + str(k) + '.txt', x_data[k], delimiter = ',')
    np.savetxt('y' + str(k) + '.txt', offset_free_data[k], delimiter = ',')




