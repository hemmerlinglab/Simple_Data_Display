import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import sys
from dunham_fit_functions import fit_dunham_coefficients, get_U_init

sys.path.append("../Analysis_Scripts/")
from data_functions import *
from aux_spectrum_functions import get_spectrum
from constants import *

###############################################
# Change global font size

font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)
###############################################

def get_spectra(folder):

    # read in data file
    arr = pickle.load( open( folder + '/data.pickle', "rb" ) )

    x_arr = arr['f_arr']
    y_arr = arr['sig_arr']
    line_info = arr['data']

    # apply wavemeter offset
    for k in range(len(x_arr)):
        x_arr[k] += line_info[k]['frequency_offset']

    # combine all data for a full spectrum
    x_spectrum = []
    y_spectrum = []
    for k in range(len(x_arr)):
        x_spectrum.extend(x_arr[k])
        y_spectrum.extend(y_arr[k])

    spectrum = np.array([x_spectrum, y_spectrum])

    spectrum = datasort(spectrum)

    result = { 'x_arr' : x_arr, 'y_arr' : y_arr, 'line_info' : line_info, 'spectrum' : spectrum }

    return result


def fit_each_line(spec, plot = False):

    x_arr = spec['x_arr']
    y_arr = spec['y_arr']
    line_info = spec['line_info']

    height_arr = []
    data = []
    line_type = []
    q_result = []

    if plot:
        #plt.figure()
        fig, (ax1, ax2) = plt.subplots(2,1)
    
    for k in range(len(x_arr)):

        vg = line_info[k]['vg']
        ve = line_info[k]['ve']
        Jg = line_info[k]['Jg']
        Je = line_info[k]['Je']

        x = x_arr[k]
        y = y_arr[k]

        # check if the line should be skipped for fitting
        if not 'skip' in line_info[k].keys():

            # fit P or R lines
            (result, xf, yf, y_residuals) = fit_gauss(x, y)

            if plot:
                ax1.plot(xf, yf, '-')
                ax1.plot(x, y, '.-')

                ax1.text(np.mean(xf), np.max(yf), str(k))

                # plot residuals
                ax2.plot(x, y_residuals)


            # create data matrix to fit Dunham coefficients
            hlp = [vg, Jg, ve, Je, result.params['p2'].value, line_info[k]['iso']]

            if Jg == 1 and Je == 1:
                # this case is a Q line but with a single Gaussian fit
                q_result.append(hlp)
            else:
                line_type.append(get_line_type(Jg, Je))
                height_arr.append(result.params['p1'].value)
                data.append(hlp)

        else:
            if Jg == 1 and Je == 1:
            # fit the Q transitions with multiple Gaussians 

                (result_multi, xf, yf, y_residuals) = fit_multi_lines(x, y)
                
                if plot:
                    ax1.plot(xf, yf, '-')
                    ax1.plot(x, y, '.-')

                    ax1.text(np.mean(xf), np.max(yf), str(k))
            
                    ax2.plot(x, y_residuals)
                
                # adding the fit result manually
                q_result = [[vg, Jg, ve, Je, result_multi.params['p2'].value, 35],\
                            [vg, Jg, ve, Je, result_multi.params['p2'].value + result_multi.params['isotope_shift'].value, 37]]
   
    fit_results = { 'line_height' : height_arr, 'line_data' : data, 'line_type' : line_type }

    if plot:
        fig.tight_layout()

    return (fit_results, q_result)



def get_rot_energy(J):

    Y01 = 0.24393004504986757
    Y02 = -2.501383745168216e-07

    Beff = h_planck * 100 * c * (Y01 + Y02 * J * (J + 1))

    eng = Beff * J * (J + 1)

    return eng


# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, my_lines, plot_fit = False):
    
    T = params['T']
    a = params['a']

    if plot_fit == False:
        x_fit = x
    else:
        no_points = 1000
        x_plot = np.linspace(np.min(x), np.max(x), no_points)
        x_fit = x_plot

    Ej = get_rot_energy(x_fit)

    Jarr = np.arange(0, 200, 1)
    Ntotal= np.sum( (2*Jarr+1) * np.exp(-get_rot_energy(Jarr)/(kB*T)) )
    
    model = a * (2*x_fit + 1)/Ntotal * np.exp(-Ej/(kB*T))

    if plot_fit == False:
        return model - data
    else:
        return (x_plot, model)



def fit_rotational_temp(x, y, my_lines = None):

    # fits a sum of gaussians to a data set
    # my_lines is a list of frequency offsets

    params = Parameters()
    
    params.add('T', value = 10.0, min = 0.0, max = 100.0, vary = True)
    params.add('a', value = 1.0, min = 0.0, max = 100.0, vary = True)
    
    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(x, y, my_lines))
    result = minner.minimize()
    
    # Store the Confidence data from the fit
    con_report = lmfit.fit_report(result.params)
    
    (x_plot, model) = fcn2min(result.params, x, y, my_lines, plot_fit = True)
    
    # get residuals
    (residuals) = fcn2min(result.params, x, y, my_lines)
    
    #:print(result.params)
    
    return (x_plot, model, result, residuals)


######################################################################################

# get data
rot_data = get_spectra('rotational_spectrum')

(fit_rot, fit_q00) = fit_each_line(rot_data, plot = True)


line_data = fit_rot['line_data']

inds = []

for k in range(len(line_data)):

    if line_data[k][-1] == 35 and line_data[k][0] == 0 and line_data[k][1] < line_data[k][3]:
        print(line_data[k])

        inds.append(k)

xs = []
ys = []
Jg = []

for k in inds:
    Jg.append(line_data[k][1])
    xs.append(line_data[k][4])
    ys.append(fit_rot['line_height'][k])

xs = np.array(xs)
ys = np.array(ys)
Jg = np.array(Jg)


#Ej = h_planck * 100 * c * 0.24 * Jg * (Jg+1)
#T = 10
#ys = (2*Jg + 1) * np.exp(-Ej/(kB*T))

(xf, yf, result, residuals) = fit_rotational_temp(Jg, ys)

print(result.params)


plt.figure()

plt.plot(Jg, ys/np.max(ys), 'o')
plt.plot(xf, yf/np.max(ys))

plt.text(2.0,0.5, "T = {0:0.2f}K".format(result.params['T'].value))

plt.xlabel('Rotational Quantum Number J')

plt.show()




