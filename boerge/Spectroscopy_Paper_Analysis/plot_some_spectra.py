import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys

sys.path.append("../Analysis_Scripts/")
from data_functions import *
from aux_spectrum_functions import get_spectrum
from energy_functions import get_scaled_dunham
from constants import *


def load_dunham(filename = 'dunham_matrices_fit.pickle'):
    
    arr = pickle.load( open( filename, "rb" ) )

    return (arr['Ug'], arr['Ue'])


def plot_spectrum(vg = 0, ve = 0, cnt_freq = 0.0):

    (Ug, Ue) = load_dunham()
    
    (Yg35, Ye35, Yg37, Ye37) = get_scaled_dunham(Ug, Ue)
    
    (nus, g35, g37) = get_spectrum(Ug, Ue, vg = vg, ve = ve, Jmax = 15, T = 10, df = 150e9)
    
    spec_35 = g35['P'] + g35['Q'] + g35['R']
    spec_37 = g37['P'] + g37['Q'] + g37['R']

    cnt_freq = np.mean(nus)
    nus = (nus - cnt_freq)/1e9
   
    a = 1.0
    
    plt.plot(nus, (spec_35 + spec_37)*a, 'k')
    



#################################################################
# Main
#################################################################


plt.figure()

plt.subplot(3,1,1)

plot_spectrum(vg = 0, ve = 0)

plt.subplot(3,1,2)

plot_spectrum(vg = 1, ve = 1)

plt.subplot(3,1,3)

plot_spectrum(vg = 2, ve = 2)

T = 100
n = np.arange(0,250)
pop = np.exp( -2*np.pi*hbar * 100*c*480 * (n+0*0.5)/(kB * T) )

print(pop/np.sum(pop))

plt.show()


