from helper import *
import numpy as np
import matplotlib.pyplot as plt



(Ug, Ue) = get_reduced_dunham()

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)

ve = 0
vg = 0
cnt_freq = 1146.330906e12

(nus, [P35, Q35, R35], [f_P, f_Q, f_R]) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, Jmax = 10, T = 10, real_amplitudes = True, isotope_abundance = 0.76)
(nus, [P37, Q37, R37], [f_P, f_Q, f_R]) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, Jmax = 10, T = 10, real_amplitudes = True, isotope_abundance = 0.24)

#print(energy(Yg35, 0, 0))
#print(energy(Ye35, 0, 0))

nus = nus - cnt_freq
nus = nus/1e9

plt.figure()
plt.plot(nus, P35, 'r')
plt.plot(nus, Q35, 'k')
plt.plot(nus, R35, 'b')

plt.plot(nus, P37, 'r--')
plt.plot(nus, Q37, 'k--')
plt.plot(nus, R37, 'b--')


plt.show()



