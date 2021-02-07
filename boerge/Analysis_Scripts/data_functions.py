import numpy as np
from configparser import ConfigParser
from os import path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import lmfit
from lmfit import Minimizer, Parameters, report_fit
from prettytable import PrettyTable
import pickle

import sys
sys.path.append("Analysis_Scripts/")
from constants import *
from energy_functions import get_scaled_dunham, get_energy


cnt_freq_00 = 1146.330000e12
cnt_freq_11 = 1145.330000e12 - 100e9


###########################################################################################
# Print complete Dunham fit report
###########################################################################################

def save_latex_table(M, title, caption = ''):

    no_col = 0
    for x in range(len(M)):
        no_col = np.max([no_col, len(M[x])])

    s = '\\begin{table}\n\\begin{tabular}{c' + 'r'*no_col + '}\n'

    s += '\\toprule\n' + title + '\\\\ \midrule \n'

    for x in range(len(M)):
        for y in range(len(M[x])):
            
            val = np.float(M[x][y])

            if np.abs(val) > 1.0:
                s += "& {0:2.3f}".format(M[x][y]) + ' '
            else:
                s += "& {0:2.3E}".format(M[x][y]) + ' '
            
        s += '\\\\\n'

    s += '\\bottomrule\n\\end{tabular}\caption{' + caption + '}\\end{table}\n'

    return s

def plot_errors(comparison, my_max = 200):

    d = {'0' : {'35' : { 'm' : [], 'p' : [] }, '37' : { 'm' : [], 'p' : [] }}, '1' : {'35' : { 'm' : [], 'p' : [] }, '37' : { 'm' : [], 'p' : [] }}}


    for k in range(1, len(comparison)):
        
        fm = np.float(comparison[k][4])
        fp = np.float(comparison[k][5])

        nu = str(comparison[k][0])
        
        iso = str(comparison[k][7])

        d[nu][iso]['m'].append(fm)
        d[nu][iso]['p'].append(fp)

    for k1 in d.keys():
        for k2 in d[k1].keys():
            for k3 in d[k1][k2].keys():
                d[k1][k2][k3] = np.array(d[k1][k2][k3])

    sub = lambda k1, k2 : d[k1][k2]['m'] - d[k1][k2]['p']

    plt.figure()

    plt.subplot(2,1,1)
    plt.plot((d['0']['35']['m']*1e12 - cnt_freq_00)/1e9, sub('0', '35')*1e6, 'ob', label = '35')
    plt.plot((d['0']['37']['m']*1e12 - cnt_freq_00)/1e9, sub('0', '37')*1e6, 'or', label = '37')

    plt.xlabel("Meas. Frequency Detuning (GHz) + {0:2.6f} THz".format(cnt_freq_00/1e12))
    plt.ylabel("Meas. - Cal. Freq. (MHz)")

    if not my_max == -1:
        plt.ylim(-my_max, my_max)

    plt.legend()

    plt.tight_layout()

    plt.subplot(2,1,2)
    plt.plot((d['1']['35']['m']*1e12 - cnt_freq_11)/1e9, sub('1', '35')*1e6, 'ob', label = '35')
    plt.plot((d['1']['37']['m']*1e12 - cnt_freq_11)/1e9, sub('1', '37')*1e6, 'or', label = '37')
    
    plt.xlabel("Meas. Frequency Detuning (GHz) + {0:2.6f} THz".format(cnt_freq_11/1e12))
    plt.ylabel("Meas. - Cal. Freq. (MHz)")
    
    if not my_max == -1:
        plt.ylim(-my_max, my_max)
    
    plt.legend()
    
    plt.tight_layout()
    
   
    return


def compare_exp_theory(data, Ug, Ue, Dg, De):
    
    (Yg35, Ye35, Yg37, Ye37) = get_scaled_dunham(Ug, Ue, Dg, De)
    
    all_result = [['vg', 'Jg', 've', 'Je', 'measured', 'predicted', 'difference', 'isotope', 'line_type']]
            
    avg_err = 0.0        
             
    for k in range(len(data)):
        
        # d = [ [gs_v, gs_J, ex_v, ex_J, freq, 35/37] ]

        vg = data[k][0]
        Jg = data[k][1]
        ve = data[k][2]
        Je = data[k][3]
        meas_freq = data[k][4]
        isotope = data[k][5]
        
        if isotope == 35:
            eng = get_energy(Yg35, Ye35, vg, Jg, ve, Je)
        elif isotope == 37:
            eng = get_energy(Yg37, Ye37, vg, Jg, ve, Je)
        else:
            print('Error')
            asd

        my_type = get_line_type(Jg, Je) 
        
        hlp = [vg, Jg, ve, Je, meas_freq/1e12, eng/1e12, (meas_freq - eng)/1e6, isotope, my_type]

        avg_err += np.abs((meas_freq - eng)/1e6)

        all_result.append(hlp)

    avg_err /= len(data)

    return (all_result, avg_err)


def make_report(data, Ug, Ue, Dg, De, latex_dunham_file = None, latex_prediction_file = None): #, vmax = 1, Jmax = 1, save_filename = 'report.txt'):
    
    (Yg35, Ye35, Yg37, Ye37) = get_scaled_dunham(Ug, Ue, Dg, De)
    
    # save all Dunham matrices
    save_matrices({ 'Ug' : Ug, 'Dg' : Dg, 'Ue' : Ue, 'De' : De, 'Yg35' : Yg35, 'Ye35' : Ye35, 'Yg37' : Yg37, 'Ye37' : Ye37 }, 'dunham_matrices_fit')


    # Latex output matrices

    if not latex_dunham_file is None:
       f = open(latex_dunham_file, 'w')
       
       s = save_latex_table(Ue, '$U (A^1\Pi)$', caption = 'Mass-reduced Dunham A-state. Units (1/cm).')
       s += save_latex_table(De, '$\Delta_\\textrm{Cl} (A^1\Pi)$', caption = 'Born-Oppenheimer Corrections A-State, Cl. Unit (1/cm).') 
       s += save_latex_table(Ye35, '$Y_{35} (A^1\Pi)$', caption = 'Dunham A-state, Cl-35. Unit (1/cm).') 
       s += save_latex_table(Ye37, '$Y_{37} (A^1\Pi)$', caption = 'Dunham A-state, Cl-37. Unit (1/cm).') 
       
       
       
       s += save_latex_table(Ug, '$U (X^1\Sigma)$', caption = 'Mass-reduced Dunham X-state. Unit (1/cm).') 
       s += save_latex_table(Dg, '$\Delta_\\textrm{Cl} (X^1\Sigma)$', caption = 'Born-Oppenheimer Corrections X-State, Cl. Unit (1/cm).') 
       s += save_latex_table(Yg35, '$Y_{35} (X^1\Sigma)$', caption = 'Dunham X-state, Cl-35. Unit (1/cm).') 
       s += save_latex_table(Yg37, '$Y_{37} (X^1\Sigma)$', caption = 'Dunham X-state, Cl-37. Unit (1/cm).') 
       
       f.write(s)
       
       f.close()

    print()
    print("*"*80)
    print("Report")
    print("*"*80)

    ## resulting dunham coefficients
    #print("\n" + "*"*40)
    #print("X-state")
    #print("*"*40 + "\n")

    #print_matrix(Ug, txt = 'Ug')
    #print_matrix(Yg35, txt = 'Yg-35')
    #print_matrix(Yg37, txt = 'Yg-37')

    #print()

    ## resulting dunham coefficients
    #print("\n" + "*"*40)
    #print("A-state")
    #print("*"*40 + "\n")

    #print_matrix(Ue, txt = 'Ue')
    #print_matrix(Ye35, txt = 'Ye-35')
    #print_matrix(Ye37, txt = 'Ye-37')

    #print()


    # table to compare measured and calculated lines 
    tab_title = ['vg', 'Jg', 've', 'Je', 'Exp. Freq. (THz)', 'Cal. Freq. (THz)', 'Diff. (MHz)', 'AlCl', 'Type']
    t = PrettyTable(tab_title)

    s = '\\begin{table}\n\\begin{tabu}{' + 'r'*(len(tab_title)) + '}\n'
    s += '\\toprule\n' # + title + '\\\\ \midrule \n'
    s += ' & '.join(tab_title) + '\\\\ \\midrule \n'

    t.float_format['Exp. Freq. (THz)'] = ".6"
    t.float_format['Cal. Freq. (THz)'] = ".6"
    t.float_format['Diff. (MHz)'] = "8.2"

    # compare experiment and theory
    all_result, avg_err = compare_exp_theory(data, Ug, Ue, Dg, De)

    for k in range(len(data)):
        
        hlp = all_result[k+1]
        t.add_row(hlp)

        # add to latex string
        hlp[4] = "{0:2.6f}".format(hlp[4])
        hlp[5] = "{0:2.6f}".format(hlp[5])
        hlp[6] = "{0:0.1f}".format(hlp[6])

        my_type = hlp[-1]

        if my_type == 'R':
            s += "\\rowfont{\\color{red}}" + " & ".join(list(map(str, hlp))) + '\\\\\n'
        elif my_type == 'P':
            s += "\\rowfont{\\color{blue}}" + " & ".join(list(map(str, hlp))) + '\\\\\n'
        else:
            s += " & ".join(list(map(str, hlp))) + '\\\\\n'

    t.add_row(['-']*6 + ['--------'] + ['-']*2)
    t.add_row(['','','','','','','avg:' + "{0:2.2f}".format(avg_err),'',''])
    
    print(t)

    s += '\\midrule\n'
    s += ' & '* 6 + ' avg: {0:2.1f}\\\\\n'.format(avg_err)
    s += '\\bottomrule\n\\end{tabu}\caption{Comparison of measured and predicted (calculated) transition frequencies.}\\end{table}\n'

    if not latex_prediction_file is None:
        f = open(latex_prediction_file, 'w')
        f.write(s)
        f.close()

    return all_result

###############################################
# Save Dunham matrix to file
###############################################

def save_matrices(Ms, filename):

    # Save matrices in one file
    f = open(filename + '.txt', 'w')

    for k in Ms.keys():
        M = Ms[k]
        f.write(k + " = [\n")
        for k in range(len(M)):
            f.write("[ ")
            f.write(",".join(map(str, M[k])))
            f.write(" ],\n")
        f.write("]\n\n")

    f.close()


    # Save each matrices in pickle file

    with open(filename + '.pickle', 'wb') as f:
        pickle.dump(Ms, f)

    return

###############################################
# Save line centers in csv file
###############################################

def spectrum2csv(data_original, save_filename = 'data_lines.txt', add_q = False, q_lines = []):

    # table to compare measured and calculated lines

    
    #t  = ",".join(['vg', 'Jg', 've', 'Je', 'Measured Fitted Freq. (THz)', 'Calculated Freq. (THz)', 'Difference (MHz)', 'AlCl', 'Type'])
    t  = ",".join(['vg', 'Jg', 've', 'Je', 'Measured Fitted Freq. (THz)', 'AlCl', 'Type'])
    t += "\n"

    if add_q:
        data = data_original.copy()
        data.append(q_lines[0])
        data.append(q_lines[1])
    else:
        data = data_original

    #print(data)

    #print(data)
    for k in range(len(data)):
        
        # d = [ [gs_v, gs_J, ex_v, ex_J, freq, 35/37] ]

        vg = data[k][0]
        Jg = data[k][1]
        ve = data[k][2]
        Je = data[k][3]
        meas_freq = data[k][4]
        isotope = data[k][5]

        my_type = get_line_type(Jg, Je) 
        
        #t += "{0}, {1}, {2}, {3}, {4:6.6f}, {5:6.6f}, {6:6.2f}, {7}, {8}\n".format(vg, Jg, ve, Je, meas_freq/1e12, eng/1e12, (eng - meas_freq)/1e6, isotope, my_type)
        t += "{0}, {1}, {2}, {3}, {4:6.6f}, {5}, {6}\n".format(vg, Jg, ve, Je, meas_freq/1e12, isotope, my_type)


    f = open(save_filename, 'w')
    f.write(t)
    f.close()
    
    #print(t)

###############################################

def print_params(params):

    print()
    for k in sorted(params.keys()):

        par = params[k]

        print("{0} : {2:12.6f} <= {1:12.6f} <= {3:12.6f}".format(k, par.value, par.min, par.max))

###############################################

def print_matrix(U, txt = None):

    if not txt is None:
        print(txt + ' = [')
    else:
        print('[')
    for k in range(len(U)):
        print(U[k], end = '')
        #for l in range(len(U[k])):
        #    print(U[k][l])
        if k < len(U):
            print(',')

    print(']\n')


###############################################

def print_dunham(U):

    print()
    print('Dunham coefficients')
    print('-'*30)
    for v in range(len(U)):
        for l in range(len(U[v])):

            val = U[v][l]

            if np.abs(val) > 1.0:
                print("U_{0}{1} = {2:15.6f}".format(v, l, val))
            else:
                print("U_{0}{1} = {2:15.6e}".format(v, l, val))

    print()
    return



###############################################
# Fit Multiple Gaussians for Q transitions
###############################################

def fcn2min_multi(params, x, data, func = None, plot_fit = False, no_of_points = 500):

    y_offset = params['p0']
    x0 = params['p2']
    w = params['p3']

    if plot_fit == False:
        x_fit = x
    else:
        x_fit = np.linspace(np.min(x), np.max(x), no_of_points)

    model = y_offset

    for k in range(2):
        model += params['a35'] * params['p_a' + str(k)] * np.exp( -(x_fit - x0 - k*params['shift_' + str(k)])**2/(2.0*w**2) )
        model += params['a37'] * params['p_a' + str(k)] * np.exp( -(x_fit - x0 - params['isotope_shift'] - k*params['shift_' + str(k)])**2/(2.0*w**2) )
    
    if plot_fit == False:
        return model - data
    else:
        return (x_fit, model)


def fit_multi_lines(x, y):

    y_min = np.min(y)
    y_max = np.max(y)
    x_mean = np.mean(x)
    x_min = np.min(x)
    x_max = np.max(x)

    params = Parameters()

    params.add('p0', value=y_min, min=-1.0, max=1.0, vary = True)
    params.add('p2', value=x_min+1.75e9, min=x_mean-10e9, max=x_mean+10e9, vary = True)
    params.add('p3', value=(x_max-x_min)/50.0, min=10e6, max=(x_max-x_min), vary = True)
    
    params.add('a35', value=0.76, min = 0.0, max = 1.0, vary = True)
    params.add('a37', value=0.24, min = 0.0, max = 1.0, vary = True)
    
    params.add('isotope_shift', value=(7.2e9 - 1.0e9), min=100e6, max=(x_max-x_min), vary = True)
    
    for k in range(2):
        params.add('p_a' + str(k), value=0.006, min=0.0, max=0.1, vary = True)
        params.add('shift_' + str(k), value=750e6, min=0e6, max=(x_max-x_min), vary = True)

    # do fit, here with leastsq model
    minner = Minimizer(fcn2min_multi, params, fcn_args=(x, y, None))
    result = minner.minimize()
    
    (xf, yf) = fcn2min_multi(result.params, x, y, plot_fit = True)
    
    (xf_res, yf_res) = fcn2min_multi(result.params, x, y, plot_fit = True, no_of_points = len(x))
    
    residuals = yf_res - y

    return (result, xf, yf, residuals)


####################################
# Fit Simple Gaussians
####################################

def fcn2min_gauss(params, x, data, func = None, plot_fit = False, no_of_points = 500):

    y_offset = params['p0']
    a = params['p1']
    x0 = params['p2']
    w = params['p3']

    if plot_fit == False:
        x_fit = x
    else:
        x_fit = np.linspace(np.min(x), np.max(x), no_of_points)

    model = y_offset + a * np.exp( -(x_fit - x0)**2/(2.0*w**2) )
    
    if plot_fit == False:
        return model - data
    else:
        return (x_fit, model)


def fit_gauss(x, y):

    y_min = np.min(y)
    y_max = np.max(y)
    x_mean = np.mean(x)
    x_min = np.min(x)
    x_max = np.max(x)

    params = Parameters()

    params.add('p0', value=y_min, min=-1.0, max=1.0, vary = True)
    params.add('p1', value=(y_max - y_min), min=0.0, max=3.0, vary = True)
    params.add('p2', value=x_mean, min=x_mean-10e9, max=x_mean+10e9, vary = True)
    params.add('p3', value=(x_max-x_min)/10.0, min=100e6, max=(x_max-x_min), vary = True)

    # do fit, here with leastsq model
    minner = Minimizer(fcn2min_gauss, params, fcn_args=(x, y, None))
    result = minner.minimize()
    
    (xf, yf) = fcn2min_gauss(result.params, x, y, plot_fit = True)
    
    (xf_res, yf_res) = fcn2min_gauss(result.params, x, y, plot_fit = True, no_of_points = len(x))
    
    residuals = yf_res - y

    return (result, xf, yf, residuals)


####################################

def get_line_type(Jg, Je):
    
    if Je == Jg:
        my_type = 'Q'
    elif Je == Jg - 1:
        my_type = 'P'
    elif Je == Jg + 1:
        my_type = 'R'

    return my_type

def get_line_color(Jg, Je):
    
    color_scheme = {'P' : 'b', 'Q' : 'k', 'R' : 'r'}
    
    return color_scheme[get_line_type(Jg, Je)]

####################################
# sort 2D data
####################################

def datasort(arr):
    
    # for 1D array
    if len(arr.shape)==1:
        return np.sort(arr)

    # for 2D array
    if len(arr.shape)==2:     

        ind = np.argsort(arr[0, :])
        arr[0, :] = arr[0, ind]
        arr[1, :] = arr[1, ind]

        return arr


####################################
# averages data
####################################

def av(arr, no_of_avg):
    
    # for 1D array
    if len(arr.shape)==1:
        hlp = np.zeros([int(arr.shape[0]/no_of_avg)])

        for k in range(len(hlp)):
            for m in range(no_of_avg):
                hlp[k] += arr[no_of_avg*k + m]

        return hlp/no_of_avg

    # for 2D array
    if len(arr.shape)==2:     

        hlp = np.zeros([int(arr.shape[0]/no_of_avg), arr.shape[1]])

        for k in range(len(hlp)):
            for m in range(no_of_avg):
                hlp[k] += arr[no_of_avg*k + m, :]

        return hlp/no_of_avg


####################################
# reads in config file
####################################

def read_in_config(f):
    
    config = ConfigParser()
    config.read(f)

    sensor_ids = config.sections()
    # make dictionary out of config

    sensors = {}

    for s in sensor_ids:
        opts = config.options(s)
        
        sensors[s] = {}
        for o in opts:
            sensors[s][o] = config.get(s, o)

    return sensors


####################################
# returns moving average
####################################

def moving_average(a, n = 3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


####################################
# combines data sets
####################################

def combine_data(x_arr, y_arr, arr = [], sort = True):

    # concatenates data arrays

    if len(arr) == 0:
        arr = range(len(x_arr))

    x = []
    y = []
    for n in arr:

        x.extend(x_arr[n])
        y.extend(y_arr[n])

    x = np.array(x)
    y = np.array(y)

    if sort:

        ind = np.argsort(x)

        x = x[ind]
        y = y[ind]

    return (x, y)


def check_which_data_folder():

    folder1 = '/home/molecules/software/data/'
    folder2 = '/Users/boerge/Software/offline_data/'
    folder3 = '/Users/johnr/Desktop/'

    if path.exists(folder1):
        datafolder = folder1
    elif path.exists(folder2):
        datafolder = folder2
    elif path.exists(folder3):
        datafolder = folder3
    else:
        print('Data path not found')
        asd

    return datafolder


def get_filename(my_date, my_time):

    datafolder = check_which_data_folder()

    basefolder = str(my_date)

    basefilename = datafolder + basefolder + '/' + basefolder + '_' + str(my_time)

    # read in config file
    config_file = basefilename + '_conf'
    conf = read_in_config(config_file)

    return basefilename, conf

##############################################
# Read in target scan image data
##############################################

def integrate_absorption(img, inter_x, inter_y, t_start, t_stop):

    abs_img = np.zeros([len(inter_x), len(inter_y)])

    for nx in range(len(inter_x)):
        for ny in range(len(inter_y)):
    
            lin_ind = nx * len(inter_y) + ny
            
            absorption = np.mean(img[lin_ind, t_start:t_stop])
            
            abs_img[nx, ny] = np.abs(absorption)

    return np.transpose(abs_img)


def integrate_patches(img, inter_x, inter_y, t, patches):

    abs_traces = []

    for k in range(len(patches)):
        abs_traces.append(np.zeros(len(t)))

    for k in range(len(patches)):

        xmin = patches[k][0][0] - patches[k][1][0]
        xmax = patches[k][0][0] + patches[k][1][0]
        ymin = patches[k][0][1] - patches[k][1][1] 
        ymax = patches[k][0][1] + patches[k][1][1]

        cnt = 0
        for nx in range(len(inter_x)):
            for ny in range(len(inter_y)):

                if      ((xmin <= inter_x[nx]) and (inter_x[nx] <= xmax)) \
                    and ((ymin <= inter_y[ny]) and (inter_y[ny] <= ymax)):
                
                    lin_ind = nx * len(inter_y) + ny
            
                    abs_traces[k] += img[lin_ind, :]
                    cnt += 1

        abs_traces[k] /= cnt

    return abs_traces


def plot_rectangle(arr, c):

    xmin = arr[0][0] - arr[1][0]
    xmax = arr[0][0] + arr[1][0]
    ymin = arr[0][1] - arr[1][1] 
    ymax = arr[0][1] + arr[1][1]
   
    plt.plot([xmin, xmax], [ymin, ymin], c)
    plt.plot([xmin, xmax], [ymax, ymax], c)
    plt.plot([xmin, xmin], [ymin, ymax], c)
    plt.plot([xmax, xmax], [ymin, ymax], c)

def plot_target_img(d, opts):

    t = d['times']
    raw_img = d['channels'][opts['channel']]

    t1 = opts['t_integrate'][0]
    t2 = opts['t_integrate'][1]

    dt = t[1] - t[0]

    t1_cut = np.where( np.abs(t - t1) <= dt )[0][0]
    t2_cut = np.where( np.abs(t - t2) <= dt )[0][0]


    img = integrate_absorption(raw_img, d['x_pos'], d['y_pos'], t1_cut, t2_cut)

    target_yields = integrate_patches(raw_img, d['x_pos'], d['y_pos'], t, opts['img_integrate'])

    plt.figure(figsize = (10,6))
    plt.subplot(1,2,1)
    plt.pcolor(d['x_pos'], d['y_pos'], img, shading = 'auto')
    #plt.colorbar()
   
    my_color_arr = 'rkbyg'

    for k in range(len(opts['img_integrate'])):
        plot_rectangle(opts['img_integrate'][k], my_color_arr[k])

    plt.xlabel('x pos')
    plt.ylabel('y pos')

    plt.gca().invert_yaxis()
    
    ax = plt.subplot(1,2,2)

    absorption_yields = []

    for k in range(len(target_yields)):
        
        # sig == 0 : 0%
        # sig == -offset : 100%
        # sig/offset - 1

        absorption = 100 * (-target_yields[k])/(d['no_absorption_level'][k] - opts['photodiode_offset'])
        
        absorption_yields.append(np.mean(absorption[t1_cut:t2_cut]))

        plt.plot(t, absorption, color = my_color_arr[k], label = opts['img_integrate'][k][2])
   
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        
        plt.ylim(-5.0,100.5)
        plt.ylabel('Absorption (%)')

    plt.legend()
    plt.title("Scan {0}_{1}".format(d['date'], d['timestamp']))
    #plt.tight_layout()

    return (img, absorption_yields)


def get_img_data(data, options):

    basefilename, conf = get_filename(data['date'], data['time'])

    # get number of averages
    no_of_avg = int(conf['scan_count']['val'])
    print('Found ' + str(no_of_avg) + ' averages.')
   
    times = read_in_file(basefilename + '_times', 1)
    
    ch0 = read_in_file(basefilename + '_ch0_arr', no_of_avg)
    posx = read_in_file(basefilename + '_posx', 1)
    posy = read_in_file(basefilename + '_posy', 1)

    # subtracting the DC offset
    offset_levels = []
    
    if options['subtract_dc_offset']:
        for k in range(ch0.shape[0]):
            
            offset_levels.append(np.mean(ch0[k, -options['offset_avg_points']:-1]))

            ch0[k, :] = ch0[k, :] - offset_levels[-1]


    inter_x = np.unique(posx)
    inter_y = np.unique(posy)

    result = {'date' : data['date'], 'timestamp' : data['time'], 'times' : times, 'x_pos' : inter_x, 'y_pos' : inter_y, 'channels' : [ch0], 'no_absorption_level' : offset_levels }

    return result


####################################
# reads in a scan
####################################

def read_in_file(filename, no_of_avg):

    x = np.genfromtxt(filename, delimiter=",")

    x = av(x, no_of_avg)

    return x

def read_laser_scan(data, options, do_print = True, frequency_domain = 'UV'):

    # data = {
    # 'date' : 20200720,
    # 'time' : 1702329_100
    # }

    # all files that Artiq produces
    # 20200807_184357_1004_act_freqs    # readout of wavemeter frequencies in THz
    # 20200807_184357_1004_ch0_arr      # transient data
    # 20200807_184357_1004_ch1_arr
    # 20200807_184357_1004_ch2_arr
    # 20200807_184357_1004_ch3_arr
    # 20200807_184357_1004_ch4_arr
    # 20200807_184357_1004_conf         # config file
    # 20200807_184357_1004_freqs        # unique setpoint frequencies in MHz
    # 20200807_184357_1004_sequence     # sequence file
    # 20200807_184357_1004_set_points   # setpoint frequencies including averages
    # 20200807_184357_1004_times        # time array

    basefilename, conf = get_filename(data['date'], data['time'])

    if do_print:
        print()
        print("-"*100)
        print('Analyzing file ... ' + str(data['date']) + '_' + str(data['time']))
        print("-"*100)
        print("Path: " + basefilename)
        print()
    
    # read in setpoints
    no_of_avg = int(conf['scan_count']['val'])
    no_of_setpoints = int(conf['setpoint_count']['val'])
    if do_print:
        print('Found {0} averages ...'.format(no_of_avg))
        print('Found {0} setpoints ...'.format(no_of_setpoints))


    set_freqs = read_in_file(basefilename + '_set_points', no_of_avg)
    
    try:
        act_freqs = read_in_file(basefilename + '_act_freqs', no_of_avg)
    except:
        print('Actual frequency file not found ...')
        act_freqs = np.copy(set_freqs)
    
    times = read_in_file(basefilename + '_times', 1)
    
    ch0 = read_in_file(basefilename + '_ch0_arr', no_of_avg)
    ch1 = read_in_file(basefilename + '_ch1_arr', no_of_avg)
    ch2 = read_in_file(basefilename + '_ch2_arr', no_of_avg)
    ch3 = read_in_file(basefilename + '_ch3_arr', no_of_avg)


    # subtracting the DC offset
    if options['subtract_dc_offset']:
        for k in range(ch0.shape[0]):
            ch0[k, :] = ch0[k, :] - np.mean(ch0[k, -options['offset_avg_points']:-1])
            ch1[k, :] = ch1[k, :] - np.mean(ch1[k, -options['offset_avg_points']:-1])
            ch2[k, :] = ch2[k, :] - np.mean(ch2[k, -options['offset_avg_points']:-1])
            ch3[k, :] = ch3[k, :] - np.mean(ch3[k, -options['offset_avg_points']:-1])

    # the factor to multiply the wavemeter frequencies with
    # sometimes we use UV, sometimes blue, sometimes IR
    if frequency_domain == 'UV':
        frequency_factor = 3.0
    elif frequency_domain == 'blue':
        frequency_factor = 2.0
    elif frequency_domain == 'IR':
        frequency_factor = 1.0
    else:
        print('Frequency domain not defined.')
    
    # get frequency scan interval in terms of absolute frequencies

    # check which laser was scanned
    if int(conf['which_scanning_laser']['val']) == 1:
        laser_offset = frequency_factor * np.float(conf['offset_laser1']['val']) * 1e12
    elif int(conf['which_scanning_laser']['val']) == 2:
        laser_offset = frequency_factor * np.float(conf['offset_laser2']['val']) * 1e12
    
    if 'wavemeter_offset' in conf.keys():
        wavemeter_offset = frequency_factor * np.float(conf['wavemeter_offset']['val']) * 1e6
    else:
        wavemeter_offset = 0.0

    # convert to the UV and add the offset
    set_freqs = frequency_factor * set_freqs * 1e6 + laser_offset + wavemeter_offset # in Hz

    act_freqs = frequency_factor * act_freqs * 1e12

    result = {}

    # apply moving time average 
    moving_avg_no = options['moving_avg_no']
    if moving_avg_no > 0:
        
        times = moving_average(times, n = moving_avg_no)

        ch0_avg = np.zeros([ch0.shape[0], len(times)])
        ch1_avg = np.zeros([ch1.shape[0], len(times)])
        ch2_avg = np.zeros([ch2.shape[0], len(times)])
        ch3_avg = np.zeros([ch3.shape[0], len(times)])
    
        for k in range(ch0.shape[0]):
            ch0_avg[k, :] = moving_average(ch0[k, :], n = moving_avg_no)
            ch1_avg[k, :] = moving_average(ch0[k, :], n = moving_avg_no)
            ch2_avg[k, :] = moving_average(ch0[k, :], n = moving_avg_no)
            ch3_avg[k, :] = moving_average(ch0[k, :], n = moving_avg_no)

        result = {'times' : times, 'set_freqs' : set_freqs, 'act_freqs' : act_freqs, 'channels' : [ch0_avg, ch1_avg, ch2_avg, ch3_avg], 'laser_offset' : laser_offset}

    else:
        
        result = {'times' : times, 'set_freqs' : set_freqs, 'act_freqs' : act_freqs, 'channels' : [ch0, ch1, ch2, ch3], 'laser_offset' : laser_offset}

    return result



#############################################################################
# plot single scan
#############################################################################

def plot_freq_scan(d, opts):

    plot_fac = 1e9
    t_max = 20

    t = d['times']
    x = d['set_freqs']
    x_act = d['act_freqs']
    y2d = d['channels'][opts['channel']]

    cnt_freq = d['laser_offset']
    
    t1 = opts['t_integrate'][0]
    t2 = opts['t_integrate'][1]

    dt = t[1] - t[0]

    t1_cut = np.where( np.abs(t - t1) <= dt )[0][0]
    t2_cut = np.where( np.abs(t - t2) <= dt )[0][0]
    t_max_plot = np.where( np.abs(t - t_max) <= dt )[0][0]

    xmin = np.min(x - cnt_freq)/plot_fac
    xmax = np.max(x - cnt_freq)/plot_fac

    # integrate between the two times
    y = np.abs(np.mean( y2d[:, t1_cut:t2_cut], axis = 1 ))


    plt.figure(figsize = (12, 6))

    plt.subplot(3,1,1)

    plt.pcolor((x - cnt_freq)/plot_fac, t[t1_cut:t_max_plot], np.transpose(y2d[:, t1_cut:t_max_plot]))

    plt.axhline(t[t2_cut], ls = '--', color = 'r')

    plt.ylabel('Time (ms)')

    plt.xlim(xmin, xmax)

    plt.tight_layout()

    plt.subplot(3,1,2)

    plt.plot((x - cnt_freq)/plot_fac, y)

    plt.ylabel('Signal (a.u.)')

    plt.xlim(xmin, xmax)
    
    plt.tight_layout()

    plt.subplot(3,1,3)

    plt.plot((x - cnt_freq)/plot_fac, (x - x_act)/1e6, '.')

    plt.xlim(xmin, xmax)

    plt.xlabel("Frequency UV (GHz) + {0:12.6f} THz".format(d['laser_offset']/1e12))
    plt.ylabel('Set - Act Freq (MHz)')

    plt.tight_layout()

    return (x, y)


def integrate_freq(d, opts):

    t = d['times']
    x = d['set_freqs']
    y2d = d['channels'][opts['channel']]
    
    cnt_freq = d['laser_offset']
    
    x1 = opts['freq_integrate'][0]
    x2 = opts['freq_integrate'][1]
 
    t1 = opts['t_integrate'][0]
    t2 = opts['t_integrate'][1]

    dt = t[1] - t[0]
    
    t1_cut = np.where( np.abs(t - t1) <= dt )[0][0]
    t2_cut = np.where( np.abs(t - t2) <= dt )[0][0]

    dx = x[1] - x[0]

    x1_cut = np.where( np.abs(x - cnt_freq - x1) <= dx )[0][0]
    x2_cut = np.where( np.abs(x - cnt_freq - x2) <= dx )[0][0]

    y = np.mean( y2d[x1_cut:x2_cut, :], axis = 0 )

    return (t, y)

def plot_time_scan(d, opts):

    plot_fac = 1e9
    t_max = 20

    t = d['times']
    x = d['set_freqs']
    x_act = d['act_freqs']
    y2d = d['channels'][opts['channel']]

    cnt_freq = d['laser_offset']
    
    x1 = opts['freq_integrate'][0]
    x2 = opts['freq_integrate'][1]
 
    t1 = opts['t_integrate'][0]
    t2 = opts['t_integrate'][1]

    dt = t[1] - t[0]
    
    t1_cut = np.where( np.abs(t - t1) <= dt )[0][0]
    t2_cut = np.where( np.abs(t - t2) <= dt )[0][0]
    t_max_plot = np.where( np.abs(t - t_max) <= dt )[0][0]

    dx = x[1] - x[0]

    x1_cut = np.where( np.abs(x - cnt_freq - x1) <= dx )[0][0]
    x2_cut = np.where( np.abs(x - cnt_freq - x2) <= dx )[0][0]
    #x_max_plot = np.where( np.abs(x - x_max) <= dx )[0][0]

    xmin = np.min(x - cnt_freq)/plot_fac
    xmax = np.max(x - cnt_freq)/plot_fac

    # integrate between the two times
    y = np.mean( y2d[x1_cut:x2_cut, :], axis = 0 )

    plt.figure(figsize = (12, 6))

    plt.subplot(1,2,1)

    plt.pcolor((x - cnt_freq)/plot_fac, t[t1_cut:t_max_plot], np.transpose(y2d[:, t1_cut:t_max_plot]))
    
    plt.axvline((x[x1_cut] - cnt_freq)/plot_fac, ls = '--', color = 'r')
    plt.axvline((x[x2_cut] - cnt_freq)/plot_fac, ls = '--', color = 'r')
    
    plt.xlabel("Frequency UV (GHz) + {0:12.6f} THz".format(d['laser_offset']/1e12))
    
    plt.ylabel('Time (ms)')

    #plt.xlim(xmin, xmax)

    plt.tight_layout()

    plt.subplot(1,2,2)

    plt.plot(t, y)

    plt.ylim(np.min(y)*1.1, np.mean(y[-10:] + np.abs(np.min(y))/10.0))

    plt.ylabel('Signal (a.u.)')

    plt.xlabel('Time (ms)')

    plt.tight_layout()

    return (t, y)


####################################################
# Cross section functions
####################################################

def get_number_of_molecules(I, I0, lamb):

    T = 10
    delta = 0.0

    my_gamma = 2*np.pi/lamb * np.sqrt(2*kB*T/(massAl+massCl_35))

    G_D = 1/(np.sqrt(np.pi) * my_gamma) * np.exp(-delta**2/my_gamma**2) 

    sigma = 1.0/4.0 * lamb**2 * my_gamma * G_D

    z = 2e-2 # cell length

    n = -1.0/(sigma * z) * np.log(I/I0)

    return n

