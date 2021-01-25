import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys

sys.path.append("../Analysis_Scripts/")
from data_functions import *
from aux_spectrum_functions import get_spectrum
from energy_functions import get_scaled_dunham
from constants import *



#################################################################
# Calculate Franck Condon
#################################################################

from scipy.special import factorial
from scipy.special import eval_hermite


def psi(x, we, n, mu):

    alpha = np.sqrt(mu * we/hbar)

    A = np.sqrt(alpha/(np.sqrt(np.pi)*(2**n)*factorial(n)))

    return A * np.exp(-0.5 * alpha**2 * x**2) * eval_hermite(n, alpha * x)#, out=None)

def integrate_fc(we1, Be1, n1, we2, Be2, n2, mu):

    xmax = 1e-9
    steps = 100000
    x = np.linspace(-xmax, xmax, steps)

    dx = x[1] - x[0]

    delta1 = np.sqrt(hbar)/(np.sqrt(2*mu*Be1))
    delta2 = np.sqrt(hbar)/(np.sqrt(2*mu*Be2))

    fc = psi(x - delta1, 2*np.pi*we1, n1, mu) * psi(x - delta2, 2*np.pi*we2, n2, mu)
    
    return np.sum(fc * dx)**2


def calc_fc(Yg, Ye, vg, ve):

    Jg = 0
    Je = 1

    # ground state Bernath
    we1 = 0.0
    for k in range(1, len(Yg)):
        # cycle through all orders of Yg
        we1 += Yg[k][0] * (vg + 0.5)**k

    Be1 = Yg[0][1] + Yg[0][2]*(Jg+1)*Jg + Yg[1][1]*(vg+1/2)

    we2 = 0.0
    for k in range(1, len(Ye)):
        # cycle through all orders of Ye
        we2 += Ye[k][0] * (ve + 0.5)**k
    
    Be2 = Ye[0][1] #+ Yg[0][2] * (n1 + 0.5)) * (n1 + 0.5)
    
    Be2 = Ye[0][1] + Ye[0][2]*(Je+1)*Je + Ye[1][1]*(ve+1/2)

    fc = integrate_fc(100*c*we1, 100*c*Be1, vg, 100*c*we2, 100*c*Be2, ve, mu_AlCl35 * amu)

    return fc

#################################################################

def load_dunham(filename = 'dunham_matrices_fit.pickle'):
    
    arr = pickle.load( open( filename, "rb" ) )

    return (arr['Ug'], arr['Ue'])

#################################################################

def load_spectrum(filename):

    d = np.genfromtxt(filename)

    x = d[:, 0]
    y = d[:, 1]

    ind = np.argsort(x)
    x = x[ind]
    y = y[ind]

    return (x, y)

#################################################################

def plot_dunham_scan(filename, vg = 0, ve = 0, cnt_freq = 0.0):

    (x, y) = load_spectrum(filename = filename)
    
    x = (x - cnt_freq)/1e9
    y = y/np.max(y)
    
    
    (Ug, Ue) = load_dunham()
    
    
    Be_arr = Ue[0][1] + np.linspace(-.01, .01, 3)
    
    for k in range(len(Be_arr)):
    
        Ue[0][1] = Be_arr[k]
    
        (Yg35, Ye35, Yg37, Ye37) = get_scaled_dunham(Ug, Ue)
    
        fc = calc_fc(Yg35, Ye35, vg, ve)
    
        (nus, g35, g37) = get_spectrum(Ug, Ue, vg = vg, ve = ve, Jmax = 15, T = 10, df = 150e9)
    
        spec_35 = g35['P'] + g35['Q'] + g35['R']
        spec_37 = g37['P'] + g37['Q'] + g37['R']
    
        nus = (nus - cnt_freq)/1e9
    
        a = -10.0
        offset = 0.0
    
        fac = 4.0
        plt.plot(nus, (spec_35 + spec_37)*a + offset + k*fac, 'k')
        plt.plot(x, y + k*fac, 'r-', markersize = 1)
    
        plt.text(-100, 1.0 + k*fac, "fc = {0:2.2f}%".format(100.0 * fc))



#################################################################
# Main
#################################################################


plt.figure()

plt.subplot(2,1,1)

plot_dunham_scan('spectrum_00.txt', vg = 0, ve = 0, cnt_freq = 1146.330000e12)

plt.subplot(2,1,2)

plot_dunham_scan('spectrum_11.txt', vg = 1, ve = 1, cnt_freq = 1145.330000e12 - 100e9)

plt.figure()

plot_spectrum(vg = 2, ve = 2)


plt.show()


