import numpy as np
import matplotlib.pyplot as plt
from helper import *


x = np.genfromtxt('freqs.txt')
y = np.genfromtxt('signal.txt')

cnt_freq = 1146.330906e12

ind = np.argsort(x)

x = x[ind]
y = y[ind]

plt.figure()
plt.plot((x - cnt_freq)/1e9, y, '.-')

(Ug, Ue) = get_reduced_dunham()

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)


ve = 0
vg = 0

T = 10

df = np.min(x) - cnt_freq

(nus, [P35, Q35, R35], [f_P, f_Q, f_R]) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, Jmax = 10, T = T, real_amplitudes = True, isotope_abundance = 0.76, df = df)
(nus, [P37, Q37, R37], [f_P, f_Q, f_R]) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, Jmax = 10, T = T, real_amplitudes = True, isotope_abundance = 0.24, df = df)

#print(energy(Yg35, 0, 0))
#print(energy(Ye35, 0, 0))

nus = nus - cnt_freq
nus = nus/1e9

#a = -0.006/4.0 * 0.5
#o = -0.001
#
#plt.plot(nus, a * P35 + o, 'r')
#plt.plot(nus, a * Q35 + o, 'k')
#plt.plot(nus, a * R35 + o, 'b')
#
#plt.plot(nus, a * P37 + o, 'r--')
#plt.plot(nus, a * Q37 + o, 'k--')
#plt.plot(nus, a * R37 + o, 'b--')


a = 0.006/2
o = 0

plt.plot(nus, a * P35 + o, 'k-')
plt.plot(nus, a * 0.75*Q35 + o, 'k-')
plt.plot(nus, a * R35 + o, 'k-')

plt.plot(nus, a * P37 + o, 'k-')
plt.plot(nus, a * 0.75*Q37 + o, 'k-')
plt.plot(nus, a * R37 + o, 'k-')


plt.plot((x - cnt_freq)/1e9, y, 'b.-', alpha = 0.2)

plt.xlim(np.min((x - cnt_freq)/1e9), np.max((x - cnt_freq)/1e9))

plt.xlabel("Frequency (GHz) + {0:3.6f} THz".format(cnt_freq/1e12))
plt.ylabel("Absorption Signal (a.u.)")

plt.show()

