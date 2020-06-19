import numpy as np
import matplotlib.pyplot as plt
from helper import *
from prettytable import PrettyTable

def get_gaussians(x, cnts, a, T = 4):

    y = np.zeros(len(x))

    for k in range(len(cnts)):

         w = get_doppler(T, cnts[k])
         y += simple_gauss(x, cnts[k], a[k], w)

    return y



def plot_vib_state(Ug = [], Ue = [], ve = 0, vg = 0, T = 4, Jmax = 20, cnt_freq = None, df = 300e9, pressure_T = 5, no_of_points = 10000, my_label = ''):


    (Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)
   
    (lines35, ampl35) = get_line_centers(Yg35, Ye35, ve, vg, Jmax = Jmax)
    (lines37, ampl37) = get_line_centers(Yg37, Ye37, ve, vg, Jmax = Jmax)

    if cnt_freq is None: 
        cnt_freq = np.mean(lines35[1])

    nus = np.linspace(cnt_freq - df, cnt_freq + df, no_of_points)
    
    g35 = []
    g37 = []
    for k in [0,1,2]:
        g35.append(get_gaussians(nus, lines35[k], 0.76*ampl35[k], T = pressure_T*T))
        g37.append(get_gaussians(nus, lines37[k], 0.24*ampl37[k], T = pressure_T*T))

    nus_plot = (nus - cnt_freq)/1e9

    plt.plot(nus_plot, g35[0], 'r-', label = my_label + '-P')
    plt.plot(nus_plot, g35[1], 'k-', label = my_label + '-Q')
    plt.plot(nus_plot, g35[2], 'b-', label = my_label + '-R')
    
    plt.plot(nus_plot, g37[0], 'r--')
    plt.plot(nus_plot, g37[1], 'k--')
    plt.plot(nus_plot, g37[2], 'b--')
    
    plt.xlim([np.min(nus_plot), np.max(nus_plot)])
    
    plt.ylabel("Signal (a.u.)")

    plt.xlabel("Frequency (GHz) + {0:3.6f} THz".format(cnt_freq/1e12))
    plt.yticks([])

    plt.legend()

    plt.tight_layout() 

    return



def plot_lines(Ug, Ue, label = '', T = 10, cnt_freq = None):

    plt.subplot(2,1,1)
    plot_vib_state(Ug = Ug, Ue = Ue, ve = 0, vg = 0, T = T, cnt_freq = cnt_freq, df = 100e9, my_label = label)
    
    plt.subplot(2,1,2)
    plot_vib_state(Ug = Ug, Ue = Ue, ve = 1, vg = 1, T = T, cnt_freq = cnt_freq, df = 100e9, my_label = label)
    
    return


def print_line(t, Yg, Ye, vg, Jg, ve, Je, label = ''):
 
    hlp = 100*c*(energy(Ye, ve, Je) - energy(Yg, vg, Jg))
    #print("v={0}/J={4} -> v'={1}/J={5} : {2:3.6f} = {3:3.6f} (IR)\n".format(ve, vg, hlp/1e12, hlp/3/1e12, Jg, Je))
    #print("v={0}/J={4} -> v'={1}/J={5} : {2:3.6f} = {3:3.6f} (IR)\n".format(ve, vg, hlp/1e12, hlp/3/1e12, Jg, Je))
    
    t.add_row([vg, Jg, ve, Je, hlp/1e12, hlp/3/1e12, label])

def list_lines(Ug, Ue, vg, ve, Jmax = 10):

    (Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)
  
    t = PrettyTable(['vg', 'Jg', 've', 'Je', 'Frequency (THz)', 'IR Freq. (THz)', 'AlCl'])

    t.float_format['Frequency (THz)'] = ".6"
    t.float_format['IR Freq. (THz)'] = ".6"

    print("*"*40)
    print("***** AlCl *****")
    print("*"*40)
    print("")
    print("***** P-Transitions *****")

    # P lines
    for Jg in range(Jmax):
        for Je in range(Jmax):
            if (Jg == Je + 1) and (Je >= 1):
                print_line(t, Yg35, Ye35, vg, Jg, ve, Je, '35')
                print_line(t, Yg37, Ye37, vg, Jg, ve, Je, '37')
    
    print(t)
    
    t = PrettyTable(['vg', 'Jg', 've', 'Je', 'Frequency (THz)', 'IR Freq. (THz)', 'AlCl'])

    t.float_format['Frequency (THz)'] = ".6"
    t.float_format['IR Freq. (THz)'] = ".6"
    
    print("")
    print("***** Q-Transitions *****")

    # Q lines
    for Jg in range(Jmax):
        for Je in range(Jmax):
            if (Jg == Je) and (Je >= 1):
                print_line(t, Yg35, Ye35, vg, Jg, ve, Je, '35')
                print_line(t, Yg37, Ye37, vg, Jg, ve, Je, '37')
 
    print(t)

    print("")
    print("***** R-Transitions *****")

    t = PrettyTable(['vg', 'Jg', 've', 'Je', 'Frequency (THz)', 'IR Freq. (THz)', 'AlCl'])

    t.float_format['Frequency (THz)'] = ".6"
    t.float_format['IR Freq. (THz)'] = ".6"
    
    # R lines
    for Jg in range(Jmax):
        for Je in range(Jmax):
            if (Jg == Je - 1) and (Je >= 1):
                print_line(t, Yg35, Ye35, vg, Jg, ve, Je, '35')
                print_line(t, Yg37, Ye37, vg, Jg, ve, Je, '37')
    
    print(t)



(Ug, Ue) = load_dunham('../20200525/Ug_fit.txt', '../20200525/Ue_fit.txt')

list_lines(Ug, Ue, 0, 0)


#plt.figure()
#plot_lines(Ug, Ue, label = 'Our fit', cnt_freq = 1146.330906e12)
#
#(Ug, Ue) = load_dunham('Ug_ram.txt', 'Ue_ram.txt')
#
#plt.figure()
#plot_lines(Ug, Ue, label = 'Ram fit')
#
#
#
#plt.show()




