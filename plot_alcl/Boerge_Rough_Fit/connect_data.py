import numpy as np
from helper import *
import matplotlib.pyplot as plt
from fit_func import *

datafolder = '/Users/boerge/Software/offline_data/'



data = [
        {'date':20200211, 'time':144344, 'ind1':97, 'ind2':102},
        {'date':20200211, 'time':144905, 'ind1':97, 'ind2':102},
        {'date':20200211, 'time':145721, 'ind1':97, 'ind2':102},
        {'date':20200213, 'time':124050, 'ind1':97, 'ind2':102},
        #{'date':20200213, 'time':130453, 'ind1':97, 'ind2':102}, # weird, no data
        {'date':20200213, 'time':132440, 'ind1':180, 'ind2':185},
        #{'date':20200213, 'time':135425, 'ind1':180, 'ind2':185}, # line overlaps
        {'date':20200213, 'time':140938, 'ind1':180, 'ind2':185},
        {'date':20200227, 'time':111627, 'ind1':177, 'ind2':185},
        {'date':20200227, 'time':121347, 'ind1':177, 'ind2':185},
        {'date':20200309, 'time':172748, 'ind1':180, 'ind2':295},
        {'date':20200309, 'time':174813, 'ind1':180, 'ind2':295},
        {'date':20200311, 'time':130625, 'ind1':180, 'ind2':295},
        ]



f_arr = []
sig_arr = []
sig_fit_arr = []
hlp_arr = []

for n in range(len(data)):

    basefolder = str(data[n]['date'])

    basefilename = datafolder + basefolder + '/' + basefolder + '_'

    f_freqs = basefilename + str(data[n]['time']) + '_set_points'
    f_ch0 = basefilename + str(data[n]['time']) + '_ch0_arr'

    config_file = basefilename + str(data[n]['time']) + '_conf'

    conf = read_in_config(config_file)

    print('Analyzing file ... ' + f_freqs)
    
    freqs = np.genfromtxt(f_freqs, delimiter=",")
    ch0 = np.genfromtxt(f_ch0, delimiter=",")

    # get number of averages
    no_of_avg = int(len(freqs)/len(np.unique(freqs)))

    print('Found ' + str(no_of_avg) + ' averages.')

    # take the averages
    freqs = av(freqs, no_of_avg)
    ch0 = av(ch0, no_of_avg)

    # subtracting the DC offset
    offset_avg_points = 5
    for k in range(ch0.shape[0]):
        ch0[k, :] = ch0[k, :] - np.mean(ch0[k, -offset_avg_points:-1])


    # get frequency scan interval in terms of absolute frequencies
    laser_offset = np.float(conf['cooling_offet']['val'])

    freq_interval = (3.0/2.0 * freqs*1e6 + 3.0 * laser_offset*1e12)

    ind1 = data[n]['ind1']
    ind2 = data[n]['ind2']
    signal = -np.mean(ch0[:, ind1:ind2], axis = 1)
    hlp = np.mean(ch0, axis = 0)

    ## fit each frequency slice for the absorption signal
    #hlp_fit = []
    #for nn in range(ch0.shape[0]):
    #
    #    a_sig = np.mean(ch0[nn:nn+1, ind1:250], axis = 0)
    #    x_sig = np.arange(len(a_sig))
    #
    #    params = Parameters()
 
    #    params.add('p0', value=0.0, min=-1.0, max=1.0, vary = False)
    #    params.add('p1', value=-1.0, min=-3.0, max=0.0, vary = True)
    #    params.add('p2', value=0.0, min=-5, max=5, vary = True)
    #    params.add('p3', value=2.3, min=0.001, max=100.0, vary = True)
    #
    #    (xf, yf, m) = fit_func(x_sig, a_sig, params, abs_sig)

    #    area = -m.params['p1']*m.params['p3']**2

    #    hlp_fit.append(area) # save amplitude
    #       

    #    plt.plot(xf, yf)
    #    plt.plot(x_sig, a_sig, 'o')


    # apply moving average
    n_mov = 3
    signal = moving_average(signal, n = n_mov)
    freq_interval = moving_average(freq_interval, n = n_mov)
    hlp = moving_average(hlp, n = n_mov)
    #hlp_fit = moving_average(hlp_fit, n = 3)

    # append all signals to arrays
    f_arr.append(freq_interval)
    sig_arr.append(signal)
    hlp_arr.append(hlp)
    #sig_fit_arr.append(hlp_fit)


f_arr = np.array(f_arr)
sig_arr = np.array(sig_arr)
hlp_arr = np.array(hlp_arr)
#sig_fit_arr = np.array(sig_fit_arr)

#ind = np.argsort(f_arr)

#f_arr = f_arr[ind]
#sig_arr = sig_arr[ind]

plt.figure()

c = 299792458

cnt_freq = 100.0 * c * 38237.00

for n in range(len(data)):
    plt.plot((f_arr[n] - cnt_freq)/1e9, sig_arr[n], 'o-', label = data[n]['time'])

plt.legend()

#plt.figure()
#
#for n in range(len(data)):
#    plt.plot((f_arr[n] - cnt_freq)/1e9, sig_fit_arr[n], 'x-', label = data[n]['time'])
#
#plt.legend()


#plt.figure()
##plt.plot(hlp_arr.transpose())
#for n in range(len(data)):
#    plt.plot(hlp_arr[n], 'o-', label = data[n]['time'])
#
#plt.legend()
#
#plt.figure()
#
#plt.plot(np.mean(hlp_arr[:, :], axis = 0))




# save arrays

import pickle

with open('data.pickle', 'wb') as f:
    pickle.dump({'f_arr' : f_arr, 'sig_arr' : sig_arr, 'data' : data}, f)


plt.show()



