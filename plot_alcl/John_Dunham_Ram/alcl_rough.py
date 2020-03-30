import numpy as np
from configparser import ConfigParser


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


def get_E(Y,v,J):
    #print(Y)
    E = 0.0
    d = 3
    for i in range(d):
        for j in range(d):
            E += Y[i][j] * (v+0.5)**i * (J*(J+1.0))**j

    return E


def get_Dunham_mat():
    molidsX, molXs = read_in_dunham_config()
    molmatX = molXs['AlCl62X_Bernath']['matrix'][:3,:3]
    molmatA = molXs['AlCl62A_Ram_Fit']['matrix']
    return molmatX,molmatA

def getP(vx,va,n):
	YX,YA = get_Dunham_mat()
	lines = []
	for i in range(n):
		lines.append((get_E(YA,va,i)-get_E(YX,vx,i+1)))
	return lines

def getQ(vx,va,n):
	YX,YA = get_Dunham_mat()
	lines = []
	for i in range(n):
		lines.append((get_E(YA,va,i)-get_E(YX,vx,i)))
	return lines


def getR(vx,va,n):
	YX,YA = get_Dunham_mat()
	lines = []
	for i in range(n):
		lines.append((get_E(YA,va,i+1)-get_E(YX,vx,i)))
	print lines


def ftw(freq):
	return freq*1e10/299792458

def wtf(wavenum):
	return wavenum*299792458/1e10




if __name__ == '__main__':
	q00 = getQ(0,0,5)
	r00 = getR(0,0,10)
	p00 = getP(0,0,10)

	print('Q00:')
	for q in range(len(q00)):
		print('{} -> {}: {}'.format(q,q,q00[q])
	print('R00:')
	for r in range(len(r00)):
		print('{} -> {}: {}'.format(r,r+1,r00[r]))
	print('P00:')
	for p in range(len(p00)):
		print('{} -> {}: {}'.format(p,p+1,p00[p]))

