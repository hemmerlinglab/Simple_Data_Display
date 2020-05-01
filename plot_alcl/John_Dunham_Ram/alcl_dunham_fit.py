import numpy as np
from configparser import ConfigParser
import lmfit
from lmfit import Minimizer, Parameters, report_fit
import matplotlib.pyplot as plt

def read_in_dunham_config(filename = 'dunham_config.ini'):
    config = ConfigParser()
    config.read(filename)
    molecule_ids = config.sections()
    molecules = {}
    
    for mi in molecule_ids:
        molmat = np.zeros((5,5))
        opts = config.options(mi)
        molecules[mi] = {}
        for o in opts:
            newval = np.float(config.get(mi,o))
            molecules[mi][o] = newval
            mollist = list(o)
            molmat[int(mollist[1]),int(mollist[2])] = newval
        molecules[mi]['matrix'] = molmat

    return molecule_ids,molecules

def get_E(Y,v,J,mat_size):
    #print(Y)
    E = 0.0
    #d = 4
    for i in range(mat_size[0]):
        for j in range(mat_size[1]):
            E += Y[i][j] * (v+0.5)**i * (J*(J+1.0))**j

    return E


def wtf(val_in_waves):
	return val_in_waves*100*299792458*1e-12


def get_line(v1,J1,v2,J2):
	mat_size = [5,5]
	return (get_E(YA,v2,J2,mat_size) - get_E(YX,v1,J1,mat_size))


mol_ids,mols = read_in_dunham_config()
YX = mols['AlCl62X_Bernath']['matrix']
YA = mols['AlCl62A_Ram_Fit']['matrix']
lines = []
for i in range(5):
	lines.append(wtf(get_line(i,0,i,0)))
	print(i,lines[i])


