import numpy as np
import datetime
import matplotlib.pyplot as plt
from configparser import ConfigParser
from os import path


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
            'S1/2-P1/2' : 377.107463380e12,
            'S1/2' : {
                '2' : -2.563005979e9,
                '1' : 4.271676631e9
                },
            'P3/2' : {
                '3' : +193.7407e6,
                '2' : -72.9112e6,
                '1' : -229.8518e6,
                '0' : -302.0738e6
                },
            'P1/2' : {
                '2' : +305.44e6,
                '1' : -509.06e6,
                }
            }
    
    Rb85 = {
            'atom' : 'Rb85',
            'S1/2-P3/2' : 384.230406373e12,
            'S1/2-P1/2' : 377.107385690e12,
            'S1/2' : {
                '3' : -1.264888516e9,
                '2' : 1.770843922e9
                },
            'P3/2' : {
                '4' : +100.205e6,
                '3' : -20.435e6,
                '2' : -83.835e6,
                '1' : -113.208e6
                },
            'P1/2' : {
                '3' : +150.659e6,
                '2' : -210.923e6,
                }

        }
   
    
    # only keeping the lines that cycle photons
    
    my_lines = []
    my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '3'))
    my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '2', Fe_co = '3', crossover = True))
    my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '2', '1', Fe_co = '3', crossover = True))
    
    
    my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '4'))
    my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '3', Fe_co = '4', crossover = True))
    my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '3', '2', Fe_co = '4', crossover = True))

    #my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '3'))
    #my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '4', Fe_co = '3', crossover = True))
    #my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '2', Fe_co = '3', crossover = True))
    #my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '1', Fe_co = '3', crossover = True))
 
    #my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '3'))
    #my_lines.append(create_transition(Rb85, 'S1/2', 'P3/2', '2', '2'))
    #my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '1', '2'))
    #my_lines.append(create_transition(Rb87, 'S1/2', 'P3/2', '1', '1'))

    my_lines = []
    my_lines.append(create_transition(Rb87, 'S1/2', 'P1/2', '1', '2'))
    my_lines.append(create_transition(Rb87, 'S1/2', 'P1/2', '1', '1', Fe_co = '2', crossover = True))
    #my_lines.append(create_transition(Rb87, 'S1/2', 'P1/2', '1', '1'))

    my_lines.append(create_transition(Rb87, 'S1/2', 'P1/2', '2', '2'))
    my_lines.append(create_transition(Rb87, 'S1/2', 'P1/2', '2', '1', Fe_co = '2', crossover = True))
    
    my_lines.append(create_transition(Rb85, 'S1/2', 'P1/2', '2', '2'))
    my_lines.append(create_transition(Rb85, 'S1/2', 'P1/2', '2', '3'))
    my_lines.append(create_transition(Rb85, 'S1/2', 'P1/2', '2', '2', Fe_co = '3', crossover = True))

    my_lines.append(create_transition(Rb85, 'S1/2', 'P1/2', '3', '2'))
    my_lines.append(create_transition(Rb85, 'S1/2', 'P1/2', '3', '3'))
    my_lines.append(create_transition(Rb85, 'S1/2', 'P1/2', '3', '2', Fe_co = '3', crossover = True))



    return my_lines




def get_img_data(mydate, mytime):

    folder1 = '/home/molecules/software/data/'
    folder2 = '/Users/boerge/software/offline_data/'

    if path.exists(folder1):
        datafolder = folder1
    elif path.exists(folder2):
        datafolder = folder2
    else:
        print('Data path not found')
        asd

    basefolder = mydate

    time_stamp = mytime

    basefilename = datafolder + basefolder + '/' + basefolder + '_'

    # get config file 
    conf = read_in_config(basefilename + time_stamp + '_conf')

    f_posx = basefilename + time_stamp + '_posx'
    f_posy = basefilename + time_stamp + '_posy'
    f_ch1 = basefilename + time_stamp + '_ch0_arr'

    posx = np.genfromtxt(f_posx, delimiter=",")
    posy = np.genfromtxt(f_posy, delimiter=",")
    ch1 = np.genfromtxt(f_ch1, delimiter=",")

    # get number of averages
    no_of_avg = int(conf['scan_count']['val'])
    print('Found ' + str(no_of_avg) + ' averages.')
    
    # average the data sets
    ch1 = av(ch1, no_of_avg)

    # subtracting the DC offset
    offset_avg_points = 5
    for k in range(ch1.shape[0]):
        ch1[k, :] = ch1[k, :] - np.mean(ch1[k, 0:offset_avg_points])

    delay_in_for_loop = np.float(conf['step_size']['val']) * 1e-6
    no_of_time_points = ch1.shape[1]
    times = np.arange(0, no_of_time_points) * (delay_in_for_loop) / 1e-3

    inter_x = np.unique(posx)
    inter_y = np.unique(posy)

    return (inter_x, inter_y, times, ch1, conf)




