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
plt.scatter((x - cnt_freq)/1e9, y,marker='.')

(Ug, Ue) = get_reduced_dunham()

(Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)


ve = 0
vg = 0

T = 5

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

LANL_62_wn = np.array([
38236.51661,
38236.02688,
38235.53393,
38235.03777,
38234.53838,
38237.48315,
38237.96642,
38238.44321,
38238.91676,
38239.38705,
38239.85406,
38240.31778,
38237.47668,
38237.46697,
38237.45401,
38237.43781,
38237.41835
])

LANL_64_wn = np.array([
38236.87762,
38236.39943,
38235.91810,
38235.43364,
38234.94604,
38237.82145,
38238.29337,
38238.75899,
38239.22145,
38239.68074,
38240.13684,
38240.58974,
38237.81516,
38237.80572,
38237.79312,
38237.77736,
38237.75843
])

def wtf(wavenum):
	return 299792458*100*wavenum

LANL_62_f = (wtf(LANL_62_wn)-cnt_freq)*1e-9
LANL_64_f = (wtf(LANL_64_wn)-cnt_freq)*1e-9
#print(LANL_62_f)





a = -0.006/2.5
o = -0.002

plt.plot(nus, a * P35 + o, 'k-')
plt.plot(nus, a * 0.75*Q35 + o, 'k-')
plt.plot(nus, a * R35 + o, 'k-')

plt.plot(nus, a * P37 + o, 'k-')
plt.plot(nus, a * 0.75*Q37 + o, 'k-')
plt.plot(nus, a * R37 + o, 'k-')

for i in range(len(LANL_62_f)):#len(LANL_62_f)):
	plt.axvline(LANL_62_f[i],ymin=0.45,ymax=0.55,color='r')
	plt.axvline(LANL_64_f[i],ymin=0.45,ymax=0.55,color='r')


plt.plot((x - cnt_freq)/1e9, y, 'b.-')#, alpha = 0.2)

plt.xlim(np.min((x - cnt_freq)/1e9), np.max((x - cnt_freq)/1e9))

plt.xlabel("Frequency (GHz) + {0:3.6f} THz".format(cnt_freq/1e12))
plt.ylabel("Absorption Signal (a.u.)")

plt.show()

