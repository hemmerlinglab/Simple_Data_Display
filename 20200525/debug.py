import numpy as np
import matplotlib.pyplot as plt
from helper import *


def get_line_freq(vg, ve, Yg, Ye, Je = 1, Jg = 1):

    hlp = 100*c*(energy(Ye, ve, Je) - energy(Yg, vg, Jg))
    print("v={0}/J={4} -> v={1}/J={5} : {2:3.6f} = {3:3.6f} (IR)\n".format(ve, vg, hlp/1e12, hlp/3/1e12, Jg, Je))


def plot_data(x, y, cnt_freq, ve = 0, vg = 0, T = 4, Jmax = 10):

    ind = np.argsort(x)

    x = x[ind]
    y = y[ind]


    plt.plot((x - cnt_freq)/1e9, y, '.-')
    #plt.plot((x - cnt_freq)/1e9, y, 'b.-', alpha = 0.2)
    
    (Ug, Ue) = get_reduced_dunham()
    
    (Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)
    
    df = np.min(x) - cnt_freq
    
    (nus, [P35, Q35, R35], [f_P, f_Q35, f_R]) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, Jmax = Jmax, T = T, real_amplitudes = True, isotope_abundance = 0.76, df = df)
    (nus, [P37, Q37, R37], [f_P, f_Q37, f_R]) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, Jmax = Jmax, T = T, real_amplitudes = True, isotope_abundance = 0.24, df = df)
    
    
    # print some transitions
    print()
    get_line_freq(0, 0, Yg35, Ye35)
    get_line_freq(0, 0, Yg37, Ye37)
    get_line_freq(1, 1, Yg35, Ye35)
    get_line_freq(1, 1, Yg37, Ye37)
    
    nus = nus - cnt_freq
    nus = nus/1e9
    
    a = 0.006/2
    o = 0
    
    plt.plot(nus, a * P35 + o, 'r-')
    plt.plot(nus, a * Q35 + o, 'k-')
    plt.plot(nus, a * R35 + o, 'b-')
    
    plt.plot(nus, a * P37 + o, 'r--')
    plt.plot(nus, a * Q37 + o, 'k--')
    plt.plot(nus, a * R37 + o, 'b--')
    
    
    plt.xlim(np.min((x - cnt_freq)/1e9), np.max((x - cnt_freq)/1e9))
    
    plt.xlabel("Frequency (GHz) + {0:3.6f} THz".format(cnt_freq/1e12))
    plt.ylabel("Absorption Signal (a.u.)")

    plt.tight_layout() 



x00 = np.genfromtxt('freqs00.txt')
y00 = np.genfromtxt('signal00.txt')

x11 = np.genfromtxt('freqs.txt')
y11 = np.genfromtxt('signal.txt')



cnt_freq00 = 1146.330906e12

cnt_freq11 = 1145.243625e12 #3*381.747875e12


plt.figure()

plt.subplot(2,1,1)
plot_data(x00, y00, cnt_freq00, ve = 0, vg = 0, T = 10)

plt.subplot(2,1,2)
plot_data(x11, y11, cnt_freq11, ve = 1, vg = 1, T = 10, Jmax = 3)



plt.show()

