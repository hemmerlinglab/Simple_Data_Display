import numpy as np
import datetime
import matplotlib.pyplot as plt
from configparser import ConfigParser



# Functions

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


def create_transition(atom, gs, es, Fg, Fe, Fe_co = '', crossover = False):

    if not crossover:
        freq = atom[gs + '-' + es] + atom[gs][str(Fg)] + atom[es][str(Fe)]
        label = atom['atom'] + ' ' + Fg + '->' + Fe
        ls = '-'
    else:
        freq = atom[gs + '-' + es] + atom[gs][str(Fg)] + (atom[es][str(Fe)] + atom[es][str(Fe_co)])/2.0
        label = atom['atom'] + ' ' + Fg + '->' + Fe + '/' + Fe_co + ' (co)'
        ls = '--'

    return [freq, label, ls]


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



def calculate_Rb_transitions():

    Rb87 = {
            'atom' : 'Rb87',
            'S1/2-P3/2' : 384.230484468e12,
            'S1/2' : {
                '2' : -2.563005e9,
                '1' : 4.27167e9
                },
            'P3/2' : {
                '3' : +193.7407e6,
                '2' : -72.9112e6,
                '1' : -229.8518e6,
                '0' : -302.0738e6
                }
        }
    
    Rb85 = {
            'atom' : 'Rb85',
            'S1/2-P3/2' : 384.230406373e12,
            'S1/2' : {
                '3' : -1.264888516e9,
                '2' : 1.770843922e9
                },
            'P3/2' : {
                '4' : +100.205e6,
                '3' : -20.435e6,
                '2' : -83.835e6,
                '1' : -113.208e6
                }
        }
    
    
    # only keeping the lines that cycle photons
    
    my_lines = []
    my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '3'))
    #my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '2'))
    #my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '1'))
    
    #my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '1', Fe_co = '2', crossover = True))
    my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '2', Fe_co = '3', crossover = True))
    my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '1', Fe_co = '3', crossover = True))
    
    my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '4'))
    #my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '3'))
    #my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '2'))
    
    #my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '2', Fe_co = '3', crossover = True))
    my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '3', Fe_co = '4', crossover = True))
    my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '2', Fe_co = '4', crossover = True))


    return my_lines


def get_data(mydate, mytime, datafolder = '/home/molecules/software/data/'):

    my_today = datetime.datetime.today()

    basefolder = mydate

    time_stamp = mytime

    basefilename = datafolder + basefolder + '/' + basefolder + '_'


    # get config file 
    conf = read_in_config(basefilename + time_stamp + '_conf')


    f_freqs = basefilename + time_stamp + '_set_points'
    
    freqs = np.genfromtxt(f_freqs, delimiter=",")
   
    # get number of averages
    no_of_avg = int(len(freqs)/len(np.unique(freqs)))
    print('Found ' + str(no_of_avg) + ' averages.')
     
    # average the data sets
    freqs = av(freqs, no_of_avg)
    
    ch = [] 
    for k in range(3):
        # define file name
        f_ch = basefilename + time_stamp + '_ch' + str(k) + '_arr'

        # read in file
        ch.append(np.genfromtxt(f_ch, delimiter=","))
        
        # average the data sets
        ch[-1] = av(ch[-1], no_of_avg)


    # average over the times since the calibration is a CW experiment 
    ch_mean = []

    for k in range(3):
        ch_mean.append(np.mean(ch[k], axis = 1))


    # subtract the offset signal
    # NOT IMPLEMENTED YET

    # normalize
    #ch_mean -= np.min(ch_mean)
    ch_mean = ch_mean/np.max(ch_mean)
    #ch_mean = 1.0 - ch_mean

    
    # get laser offset
    laser_offset = np.float(conf['setpoint_offset']['val'])

    freqs = freqs * 1e6 + laser_offset * 1e12

    return (freqs, ch_mean)


