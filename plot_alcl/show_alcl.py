import numpy as np
import datetime
import matplotlib.pyplot as plt
import os
import glob
import sys
import fit_line as fl
from conf_reader import read_config
#import fit_dunham as fd

c = 299792458 # m/s
kb =  1.38064852e-23 # m^2 kg s^-2 K^-1 (Boltzmann)
amu = 1.66053892e-27 # kg/amu
m = (35+27)*amu # kg (AlCl)


# hlp is helper variable
def plot_sub(k, times, nus, ch, myfontsize = 8):

    plt.subplot(3,2,k)

    plt.pcolor(times, nus, ch)

    plt.xlabel('Time (ms)', fontsize = myfontsize)
    plt.ylabel('Frequencies (MHz)', fontsize = myfontsize)
    plt.tick_params(labelsize=myfontsize,direction='in')

    plt.gca().invert_yaxis()


    plt.tight_layout()

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


def do_the_thing(basefilename,time_stamp,delta1,delta2):

    time_stamp = str(time_stamp)

    conf_file = basefilename + time_stamp + '_conf'
    conf_data = read_config(conf_file)

    #print(conf_data['step_size']['val'])

    f_freqs = basefilename + time_stamp + '_set_points'
    f_ch0 = basefilename + time_stamp + '_ch0_arr'
    f_ch1 = basefilename + time_stamp + '_ch1_arr'
    f_ch2 = basefilename + time_stamp + '_ch2_arr'
    f_ch3 = basefilename + time_stamp + '_ch3_arr'
    f_ch4 = basefilename + time_stamp + '_ch4_arr'

    f_ch0s = basefilename + time_stamp + '_ch0_arr_slow'
    f_ch1s = basefilename + time_stamp + '_ch1_arr_slow'
    f_ch2s = basefilename + time_stamp + '_ch2_arr_slow'
    f_ch3s = basefilename + time_stamp + '_ch3_arr_slow'
    f_ch4s = basefilename + time_stamp + '_ch4_arr_slow'



    freqs = np.genfromtxt(f_freqs, delimiter=",")
    ch0 = np.genfromtxt(f_ch0, delimiter=",")
    ch1 = np.genfromtxt(f_ch1, delimiter=",")
    ch2 = np.genfromtxt(f_ch2, delimiter=",")
    ch3 = np.genfromtxt(f_ch3, delimiter=",")
    ch4 = np.genfromtxt(f_ch4, delimiter=",")

    ch0s = np.genfromtxt(f_ch0s, delimiter=",")
    ch1s = np.genfromtxt(f_ch1s, delimiter=",")
    ch2s = np.genfromtxt(f_ch2s, delimiter=",")
    ch3s = np.genfromtxt(f_ch3s, delimiter=",")
    ch4s = np.genfromtxt(f_ch4s, delimiter=",")



    # get number of averages
    no_of_avg = int(len(freqs)/len(np.unique(freqs)))

    print('Found ' + str(no_of_avg) + ' averages.')

    freqs = av(freqs, no_of_avg)
    ch0 = av(ch0, no_of_avg)
    ch1 = av(ch1, no_of_avg)
    ch2 = av(ch2, no_of_avg)
    ch3 = av(ch3, no_of_avg)
    ch4 = av(ch4, no_of_avg)

    ch0s = av(ch0s, no_of_avg)
    ch1s = av(ch1s, no_of_avg)
    ch2s = av(ch2s, no_of_avg)
    ch3s = av(ch3s, no_of_avg)
    ch4s = av(ch4s, no_of_avg)



    offset_freq = float(conf_data['cooling_offet']['val'])
    #offset_freq = 382.08140

    #nus = 2*freqs

    delay_in_for_loop = float(conf_data['step_size']['val'])*1e-6
    #delay_in_for_loop = 100e-6
    #delay_in_for_loop = 50e-6
    no_of_time_points = ch1.shape[1]
    times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3




    #time_cut1 = 182
    #time_cut2 = time_cut1+10

    #time_cut1 = 98
    #time_cut2 = time_cut1+10

    time_cut1 = int(8.4e-3/delay_in_for_loop) + 14 + delta1
    #print(time_cut1)
    time_cut2 = time_cut1 + delta2

    color_plot = ch0[:, time_cut1:time_cut2]

    color_plot -= np.min(np.min(color_plot))

    color_plot /= np.max(np.max(color_plot))

    ### 2D PLOTS
    # if delta1 == 0:
    #     plt.figure()
    #     plt.pcolor(times[time_cut1:time_cut2], 3*freqs, color_plot)
    #     plt.xlabel('Time (us)')
    #     plt.ylabel('Relative UV frequency (MHz)')

    #     plt.colorbar()

        

    signal = np.mean(ch0[:, time_cut1:time_cut2], axis = 1)
    signal = signal/np.max(np.abs(signal))

    (xfit, yfit, result) = fl.fit_line(freqs, signal)

    cnt_peak = result.params['x0'].value

    #print('a: ',result.params['a'].value)
    #print('w: ',result.params['w'].value)
    #print('x0: ',result.params['x0'].value)
    #print('y_offset',result.params['y_offset'].value)

    abs_cnt_peak = 3 * ((offset_freq * 1e12 + cnt_peak * 1e6)/1e12) # in THz

    abs_cnt_peak_wavenumber = abs_cnt_peak * 1e12/100.0/c


    ### FIT PLOTS
    # if delta1 == 0:
    #     plt.figure()
    #     plt.plot(3*freqs, 1 - signal)
    #     plt.plot(3*xfit, 1 - yfit, 'r-')

    #     plt.ylim(0, .8)

    #     plt.ylabel('Absorption signal (a.u.)')
    #     plt.xlabel('Frequency (MHz)')
    #     plt.title(time_stamp)

    #     plt.text(-400, 0.4, "Center peak frequency: \n\n     {0:9.6f} THz \n = {1:9.4f} 1/cm".format(abs_cnt_peak, abs_cnt_peak_wavenumber))


    wid = result.params['w'].value

    return abs_cnt_peak,wid,freqs,signal,offset_freq



if __name__ == '__main__':
    my_today = datetime.datetime.today()

    #datafolder = '/Users/boerge/skynet/molecule_computer/'
    #datafolder = '/home/molecules/software/data/'

    datafolder = '\\Users\\John\\Desktop\\'

    #datafolder = '/Users/boerge/software/data/Data/molecule_computer/'

    #datafolder = '/Users/johnr/Desktop/'

    stamps = [124050,130453,132440,134013,135425,140938]
    delta1s = np.linspace(0,10,2)
    #print(delta1s)
    delta2 = 10

    basefolder = '20200213'


    #time_stamp = str(stamps[0])


    #basefilename = datafolder + basefolder + '/' + basefolder + '_'    
    basefilename = datafolder + basefolder + '\\' + basefolder + '_'

    #do_the_thing(basefilename,time_stamp)
 
    centers = np.zeros((len(stamps),len(delta1s)))
    widths = np.zeros((len(stamps),len(delta1s)))
    all_freqs = [[] for i in range(len(delta1s))]
    all_signal = [[] for i in range(len(delta1s))]

    for i in range(len(stamps)):
        for j in range(len(delta1s)):
            print('Loading data {} from {} - {}...'.format(stamps[i],int(delta1s[j]),int(delta1s[j])+delta2),end='')
            centers[i][j] , widths[i][j] , new_freqs, new_sigs, new_off = do_the_thing(basefilename,stamps[i],int(delta1s[j]),int(delta1s[j])+delta2)
            #print(new_off*3e12+new_freqs*1e6)
            all_freqs[j] = np.concatenate((all_freqs[j],new_off*3e12+new_freqs*1e6),axis=0)
            all_signal[j] = np.concatenate((all_signal[j],new_sigs),axis=0)



    #print(widths)

    #print(np.mean(centers,axis=1))
    #print(widths[:][0])

    avg_cnts = np.mean(centers,axis=1)

    Ts = np.zeros(np.shape(widths))

    for n in range(len(stamps)):
        Ts[n] = (m*c**2*(widths[n]*1e6/(avg_cnts[n]*1e12))**2)/(8*kb*np.log(2))


    plt.figure()
    for l in range(len(stamps)):
        plt.subplot(3,2,l+1)
        plt.plot(delta1s,Ts[l])
        plt.title(avg_cnts[l])
        plt.xlabel('Time Steps')
        plt.ylabel('T (K)')

    print(all_freqs)


    plt.figure()
    plt.scatter(all_freqs[0],all_signal[0],marker='.',color='b')
    plt.scatter(all_freqs[1],all_signal[1],marker='.',color='r')
    plt.title('Full Spectrum')
    plt.xlabel('Frequency (THz)')
    plt.ylabel('Abs (arb.)')







    plt.show()