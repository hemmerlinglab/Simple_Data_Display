import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import sys

sys.path.append("../Analysis_Scripts/")
from constants import *
from data_functions import *
from aux_spectrum_functions import get_transition_energies, get_spectrum

from dunham_fit_functions import fit_dunham_coefficients, get_U_init

from sequence_fitter import do_sequence_fit

###############################################
# Change global font size

font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)
###############################################



# General flags

INCLUDE_Q_LINES = True
FIT_WHICH_ISOTOPES = 35 # can be 35, 37 or 0 (both)
SET_BOB_ZERO = True

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

###############################################

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

        hlp_isotope = line_info[k]['iso']
        hlp_line_type = get_line_type(Jg, Je)
        
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
            hlp = [vg, Jg, ve, Je, result.params['p2'].value, hlp_isotope]

            if Jg == 1 and Je == 1:
                # this case is a Q line but with a single Gaussian fit
                q_result.append(hlp)

                # uncomment these three lines to include the fitting of the Q11-lines
                if INCLUDE_Q_LINES:
                    
                    if True: #(hlp_isotope == FIT_WHICH_ISOTOPES or FIT_WHICH_ISOTOPES == 0):
                        line_type.append(get_line_type(Jg, Je))
                        height_arr.append(result.params['p1'].value)
                        data.append(hlp)

            else:
                if (hlp_isotope == FIT_WHICH_ISOTOPES or FIT_WHICH_ISOTOPES == 0):
                    line_type.append(hlp_line_type)
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


###############################################

def combine_fit_results(fit_arr):

    fit_all = { 'line_height' : [], 'line_data' : [], 'line_type' : [] }
    for k in range(len(fit_arr)):
        
        for key in fit_all.keys():
            fit_all[key].extend(fit_arr[k][key])
    
    return fit_all

###############################################

def plot_spectrum(spec, fits, style = '.-', cnt_freq = 0.0, filename = 'spectrum', my_color = None, do_moving_average = False, mov_n = 3):

    x_arr = spec['x_arr']
    y_arr = spec['y_arr']
    line_type = fits['line_type']

    line_data = fits['line_data']

    full_freqs = []
    full_signal = []

    for k in range(len(x_arr)):
    
        x = x_arr[k]
        y = y_arr[k]
    
        full_freqs.extend(list(x))
        full_signal.extend(list(y))
    
        #if line_type[k] == 'Q':
        #    if my_color is None:
        #        color = 'k'
        #    else:
        #        color = my_color[0]
        #elif line_type[k] == 'R':
        #    if my_color is None:
        #        color = 'r'
        #    else:
        #        color = my_color[1]
        #elif line_type[k] == 'P':
        #    if my_color is None:
        #        color = 'b'
        #    else:
        #        color = my_color[2]
        color = 'r'

        if do_moving_average:
            x = moving_average(x, n = mov_n) 
            y = moving_average(y, n = mov_n) 

        plt.plot((x - cnt_freq)/1e9, 2*60*y, style, color = color)

    #  mark the lines that were used in the fitting process
    for k in range(len(line_data)):
        plt.axvline((line_data[k][4] - cnt_freq)/1e9, ymin = 0, ymax = 1, ls = '--', color = 'r')

    # plot zero line
    plt.axhline(0, ls = '--', color = 'k')

    plt.ylabel("Signal (a.u.)")

    plt.xlabel("Frequency (GHz) + {0:3.6f} THz".format(cnt_freq/1e12))
    plt.yticks([])

    plt.tight_layout() 

    #plt.grid()


 
    # save spectrum 
    
    full_spectrum = []
    # save full spectrum
    for k in range(len(full_freqs)):
        hlp = [full_freqs[k], full_signal[k]]
    
        full_spectrum.append(hlp)
   
    full_spectrum = np.array(full_spectrum)
    
    
    np.savetxt(filename + '.txt', full_spectrum, delimiter = ' ')


###############################################

def plot_prediction(x, y, fac = 1.0, offset = 0.0, style = '-', cnt_freq = 0.0, color = None):

    if color is None:
        c1 = 'b'
        c2 = 'k'
        c3 = 'r'
    else:
        c1 = color[0]
        c2 = color[1]
        c3 = color[2]
    
    plt.plot((x - cnt_freq)/1e9, fac*y['P'] + offset, style, color = c1)
    plt.plot((x - cnt_freq)/1e9, fac*y['Q'] + offset, style, color = c2)
    plt.plot((x - cnt_freq)/1e9, fac*y['R'] + offset, style, color = c3)


###############################################

def plot_branches(fit_data, cnt_freq, Ug, Ue, Dg, De, vg = 0, ve = 0):

    d = fit_data['line_data']

    (Yg35, Ye35, Yg37, Ye37) = get_scaled_dunham(Ug, Ue, Dg, De)
    
    (l_P35, l_Q35, l_R35) = get_transition_energies(Yg35, Ye35, vg, ve, Jmax = 8)
    (l_P37, l_Q37, l_R37) = get_transition_energies(Yg37, Ye37, vg, ve, Jmax = 8)

    plt.figure()

    for k in range(len(d)):

        Jg = d[k][1]
        Je = d[k][3]
        freq = d[k][4]

        my_color = get_line_color(Jg, Je)

        plt.plot(Jg, (freq - cnt_freq)/1e9, my_color + 'x', markersize = 10)

    
    plt.plot(l_P35[:, 1], (l_P35[:, 4] - cnt_freq)/1e9, 'b')
    plt.plot(l_P37[:, 1], (l_P37[:, 4] - cnt_freq)/1e9, 'b--')
    #plt.plot(l_Q[:, 1], (l_Q[:, 4] - cnt_freq)/1e9, 'k')
    plt.plot(l_R35[:, 1], (l_R35[:, 4] - cnt_freq)/1e9, 'r')
    plt.plot(l_R37[:, 1], (l_R37[:, 4] - cnt_freq)/1e9, 'r--')

    plt.xlabel('JX')
    plt.ylabel('Frequency detuning (GHz)')

    return


########################################################################################################################################
# Main
########################################################################################################################################

cnt_freq_00 = 1146.330000e12
cnt_freq_11 = 1145.330000e12 - 100e9


# get data
v00 = get_spectra('v00')
v11 = get_spectra('v11')


# fit each line to a gaussian
(fit_pr00, fit_q00) = fit_each_line(v00, plot = False)
(fit_pr11, fit_q11) = fit_each_line(v11, plot = False)

# save the data
spectrum2csv(fit_pr00['line_data'], save_filename = 'data_v0v0.csv', add_q = True, q_lines = fit_q00)
spectrum2csv(fit_pr11['line_data'], save_filename = 'data_v1v1.csv', add_q = True, q_lines = fit_q11)

##################################
# Fit Dunham coefficients
##################################

# only use P and R lines to fit spectrum
fit_all = combine_fit_results([fit_pr00, fit_pr11])

# save fit_all dictionary to simulate error
with open('fitted_data.pickle', 'wb') as f:
        pickle.dump(fit_all, f)


#do_sequence_fit([fit_pr00, fit_pr11])

(Ug_init, Ue_init, Dg_init, De_init) = get_U_init(set_bob_zero = SET_BOB_ZERO)

(Ug, Ue, Dg, De) = fit_dunham_coefficients(fit_all['line_data'], Ug_init, Ue_init, Dg_init, De_init, vary_groundstate = False)


my_ext = str(FIT_WHICH_ISOTOPES)
comparison = make_report(fit_all['line_data'], Ug, Ue, Dg, De, latex_dunham_file = 'latex_dunham' + my_ext, latex_prediction_file = 'latex_prediction' + my_ext)



## plot branches of P- and R-lines
#plot_branches(fit_pr00, cnt_freq_00, Ug, Ue, Dg, De)
#plot_branches(fit_pr11, cnt_freq_11, Ug, Ue, Dg, De, vg = 1, ve = 1)


##################################
# Plot errors
##################################

plot_errors(comparison, my_max = 200)

##################################
## Plot spectrum and theory
##################################

(nus_00, y35_00, y37_00) = get_spectrum(Ug, Ue, Dg, De, vg = 0, ve = 0, Jmax = 15, df = 120e9, T = 6)

(nus_11, y35_11, y37_11) = get_spectrum(Ug, Ue, Dg, De, vg = 1, ve = 1, Jmax = 15, df = 120e9, T = 6)

#nus_11 = nus_00
#y35_11 = y35_00
#y37_11 = y37_00

# moving averages for display
mov_avg = True
mov_n = 5

plt.figure()

plot_spectrum(v00, fit_pr00, style = '.-', cnt_freq = cnt_freq_00, filename = 'spectrum_00', my_color = ['r'] * 3, do_moving_average = mov_avg, mov_n = mov_n)
plot_prediction(nus_00, y35_00, fac = -2, offset = 0.0, cnt_freq = cnt_freq_00, color = ['k'] * 3)
plot_prediction(nus_00, y37_00, fac = -4, offset = 0.0, style = '--', cnt_freq = cnt_freq_00, color = ['k'] * 3)


plt.figure()

plot_spectrum(v11, fit_pr11, style = '.-', cnt_freq = cnt_freq_11, filename = 'spectrum_11', my_color = ['r'] * 3, do_moving_average = mov_avg, mov_n = mov_n)
plot_prediction(nus_11, y35_11, fac = -2, offset = 0.0, cnt_freq = cnt_freq_11, color = ['k'] * 3)
plot_prediction(nus_11, y37_11, fac = -4, offset = 0.0, style = '--', cnt_freq = cnt_freq_11, color = ['k'] * 3)

#plt.xlim(-65, 106)


plt.show()



print("Shift Cl-35: {0:0.3f} 1/cm".format( Ue[0][0] * mass_e/(massCl_35 * amu) * De[0][0] ))
print("Shift Cl-37: {0:0.3f} 1/cm".format( Ue[0][0] * mass_e/(massCl_37 * amu) * De[0][0] ))
print("Diff Cl-37: {0:0.3f} 1/cm".format( Ue[0][0] * mass_e*(1/(massCl_37 * amu) - 1/(massCl_35 * amu)) * De[0][0] ))

print()

print("Shift Cl-35: {0:0.3f} 1/cm".format( Ue[1][0] * mass_e/(massCl_35 * amu) * De[1][0] ))
print("Shift Cl-37: {0:0.3f} 1/cm".format( Ue[1][0] * mass_e/(massCl_37 * amu) * De[1][0] ))
print("Diff Cl-37: {0:0.3f} 1/cm".format( Ue[1][0] * mass_e*(1/(massCl_37 * amu) - 1/(massCl_35 * amu)) * De[1][0] ))

print()


print("Shift Cl-35: {0:0.3f} 1/cm".format( Ue[2][0] * mass_e/(massCl_35 * amu) * De[2][0] ))
print("Shift Cl-37: {0:0.3f} 1/cm".format( Ue[2][0] * mass_e/(massCl_37 * amu) * De[2][0] ))
print("Diff Cl-37: {0:0.3f} 1/cm".format( Ue[2][0] * mass_e*(1/(massCl_37 * amu) - 1/(massCl_35 * amu)) * De[2][0] ))






# to make nice circles and broken axis: 




# from brokenaxes import brokenaxes
#broken_xlims = []
#for k in range(len(plot_x)):
#    
#    avg_cnt = (plot_x[k][0] + plot_x[k][-1])/2.0
#    offset = 1.5
#    hlp = [avg_cnt - offset, avg_cnt + offset]
#    broken_xlims.append(hlp)
#
#broken_xlims.sort(key = lambda x : x[0])
#
#
#
#line_pos = np.sort(line_pos)
#
#my_color = ['r', 'k', 'b', 'g', 'y']
#
#bax = brokenaxes(xlims=broken_xlims) #, subplot_spec=sps1)
#
#for k in range(len(x_data)):
#
#    
#    bax.scatter(plot_x[k], offset_free_data[k], color = my_color[k], marker = 'o', fc = 'w', lw = 2.0)
#
#

