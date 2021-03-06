import numpy as np
from configparser import ConfigParser
import lmfit
from lmfit import Minimizer, Parameters, report_fit
import matplotlib.pyplot as plt

def read_in_ram_config(filename = 'ram_config.ini'):
    config = ConfigParser()
    config.read(filename)
    line_ids = config.sections()
    lines = {}

    for li in line_ids:
        print('Loading {}...'.format(li),end='\r')
        opts = config.options(li)
        lines[li] = {}
        for o in opts:
            lines[li][o] = config.get(li,o)

        lines[li]['data'] = read_data_file(li+'.data')
        
        lines[li]['JX'] = np.arange(len(lines[li]['data']))+int(lines[li]['firstj'])

        lilist = list(li)

        if lilist[0] == 'Q':
            lines[li]['JA'] = lines[li]['JX']
        elif lilist[0] == 'P':
            lines[li]['JA'] = []
            for num in lines[li]['JX']:
                try:
                    lines[li]['JA'].append(num - 1)
                except:
                    lines[li]['JA'].append('NAN')
            lines[li]['JA'] = np.array(lines[li]['JA'])
        else:
            print('AAAAAAAHHHHHHH!!!!!!!')

        if int(lilist[1]) == 1:
            vx = int(lilist[1])*10+int(lilist[2])
            if int(lilist[3]) == 1:
                va = 10+int(lilist[4])
            else:
                va = int(lilist[3])
        else:
            vx = int(lilist[1])
            if int(lilist[2]) == 1:
                va = 10+int(lilist[3])
            else:
                va = int(lilist[2])

        lines[li]['JA'] += 3
        lines[li]['JX'] += 3
        lines[li]['vX'] = va
        lines[li]['vA'] = vx

    return line_ids,lines


def read_data_file(filename):
    baseval = 0
    dat = []
    f = open('ram/'+filename,'r')
    fl = f.readlines()
    for x in fl:
        newval = np.float(x)
        if newval > 30000:
            dat.append(newval)
            baseval = np.floor(newval/100)*100
        elif newval == 0.0:
            dat.append('NAN')
        else:
            dat.append(baseval + newval)

    return np.array(dat)


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


def make_Dunham_mat(params,mat_size):
    molidsX, molXs = read_in_dunham_config()
    molmatX = molXs['AlCl62X_Bernath']['matrix'][:mat_size[0],:mat_size[1]]
    molmatA = np.zeros((mat_size[0],mat_size[1]))
    for p in params:
        plist = list(p)
        if plist[3] == 'X':
            pass
            #molmatX[int(plist[1]),int(plist[2])] = params[p]
        elif plist[3] == 'A':
            try:
                molmatA[int(plist[1]),int(plist[2])] = params[p]
            except:
                pass
        else:
            print('OH NO!!!')
    return molmatX,molmatA


def fcn2min(params, x, data, get_fit = False):
    mat_size = [5,5] # max_dun_indx - 1  (USE THIS LINE TO CHANGE ORDER)
    YX,YA = make_Dunham_mat(params,mat_size)
    
    model = []

    if get_fit == False:
        for i in range(len(data)):
            model.append(get_E(YA, x[0][i], x[1][i], mat_size) - get_E(YX, x[2][i], x[3][i], mat_size))

        return np.array(model-data)

    else:
        for i in range(len(data)):
            model.append(get_E(YA, x[0][i], x[1][i], mat_size) - get_E(YX, x[2][i], x[3][i], mat_size))

        return np.array(model)


def fit_dunham(q,d):
    print('Starting fit...')
    molids,mols = read_in_dunham_config()
    Ys = mols['AlCl62X_Bernath'].keys()
    allowed_Ys = ['y00','y01','y10','y11','y20','y21','y22','y12','y02']
    params = Parameters()
    for p in Ys:
        if p != 'matrix':
            #if p in allowed_Ys:
                if p == 'y00':
                    params.add(p+'A', value = 0.0, vary = True)
                elif p == 'y10':
                    params.add(p+'A', value = 0.0, max = 600.0, vary = True)
                elif p == 'y20':
                    params.add(p+'A', value = 0.0, vary = True)
                else:
                    params.add(p+'A', value = 0.0, vary = True)

            #params.add(p+'X', value = 0.0, min = -500, max = 500, vary = True)
                
    
    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(q,d))
    result = minner.minimize()
    
    # Store the Confidence data from the fit
    con_report = lmfit.fit_report(result.params)
    model = fcn2min(result.params,q,d,get_fit=True)

    return (result.params, params, con_report, model)



def compile_data():
    ids,trans = read_in_ram_config()

    all_quant = [[],[],[],[]]
    all_data = []

    for ln in ids:
        for i in range(len(trans[ln]['JX'])):
            new_dat = trans[ln]['data'][i]
            if new_dat != 'NAN':
                all_quant[0].append(trans[ln]['vA'])
                all_quant[1].append(trans[ln]['JA'][i])
                all_quant[2].append(trans[ln]['vX'])
                all_quant[3].append(trans[ln]['JX'][i])
                all_data.append(np.float(new_dat))

    return np.array(all_quant),np.array(all_data)

def update_config(result_params,dunham_config='dunham_config.ini'):
    config_final = ConfigParser()
    config_final.read(dunham_config)
    rps = result_params
    for par in rps:
        config_final.set('AlCl62A_Ram_Fit',par[:3],str(rps[par].value))

    with open('dunham_config.ini','w') as configfile:
        config_final.write(configfile)

if __name__ == '__main__':
    q,data = compile_data()
    qs = q[:,0::10]
    dats = data[0::10]

    rps, ps, cr, dplot = fit_dunham(q,data)
    update_config(rps)
    print(rps)

    plt.figure()
    plt.plot(data)
    plt.plot(dplot)
    plt.ylabel('Line Energy (wavenumber)')
    plt.xlabel('Data Index (arb)')
    plt.show()
