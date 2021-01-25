import sys
sys.path.append("../Analysis_Scripts/")

from constants import *
from data_functions import *
import matplotlib.pyplot as plt
import matplotlib

font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)



options = { 'offset_avg_points' : 5, 'moving_avg_no' : 0, 'subtract_dc_offset' : False }

plot_options = { 't_integrate' : [10.685, 14.1], 'freq_integrate' : [300e6, 600e6], 'channel' : 0 }

data = [
        {'date':20200828, 'time':151418}, # Yag on, absorption on
        {'date':20200828, 'time':152115}, # Yag on, absorption blocked
        {'date':20200828, 'time':152718} # Yag blocked, absorption on
        ]


plot_options = { 't_integrate' : [10.685, 14.1], 'freq_integrate' : [100e6, 700e6], 'channel' : 0 }

data = [
        {'date':20200518, 'time':132850}, # Yag on, absorption on
        {'date':20200518, 'time':133754}, # Yag on, absorption blocked
        ]



result = []

for k in range(len(data)):
    result.append(read_laser_scan(data[k], options))

x_arr = []
y_arr = []

for k in range(len(data)):
    #(x, y) = plot_time_scan(result[k], plot_options)
    (x, y) = integrate_freq(result[k], plot_options)
    
    x = moving_average(x, n = 5)
    y = moving_average(y, n = 5)

    x_arr.append(x)
    y_arr.append(y)

#plt.figure()
#
#
#plt.plot(x_arr[0], y_arr[0], label = 'Yag on, absorption on')
#plt.plot(x_arr[1], y_arr[1], label = 'Yag on, absorption blocked')
##plt.plot(x_arr[2], y_arr[2], label = 'Yag blocked, absorption on')
#
#plt.legend()
#
#plt.ylim(-0.005, 0.005)

t = x_arr[0]

abs_signal = 100.0 * (y_arr[0] - np.mean(y_arr[1]))/(np.mean(y_arr[0][-10:]) - np.mean(y_arr[1]))
bg_signal = 100.0 * (y_arr[1] - np.mean(y_arr[1]))/(np.mean(y_arr[0][-10:]) - np.mean(y_arr[1]))


n = get_number_of_molecules(abs_signal, np.mean(abs_signal[-10:]), lamb = c/result[0]['laser_offset'])

print(n)

plt.figure()

plt.plot(t, abs_signal)

plt.ylim(0,105)

plt.xlabel('Time (ms)')
plt.ylabel('Transmission (%)')

plt.figure()

#plt.semilogy(t, n)
plt.plot(t, n)

#plt.ylim(0,105)

plt.xlabel('Time (ms)')
plt.ylabel('Number of Molecules')



plt.show()



