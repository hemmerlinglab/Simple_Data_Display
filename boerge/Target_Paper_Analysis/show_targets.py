import sys
import matplotlib.pyplot as plt
import matplotlib

sys.path.append("../Analysis_Scripts/")
from data_functions import *

from lmfit import Minimizer, Parameters, report_fit

font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)




# define objective function: returns the array to be minimized
# function to minimize
def fcn2min(params, x, data, plot_fit = False):
    
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value

    if plot_fit == False:
        x_fit = np.array(x)
    else:
        no_points = 1000
        x_plot = np.linspace(np.min(x), np.max(x), no_points)
        x_fit = x_plot



    model = a * (1 - np.exp(-b * x_fit)) * (c - x_fit)

    if plot_fit == False:
        return model - data
    else:
        return (x_plot, model)


def fit_theory(x, y, vary = True):
    
    # fits a sum of gaussians to a data set
    # my_lines is a list of frequency offsets

    params = Parameters()
    
    #params.add('a', value = np.min(y), min = 0*np.min(y), max = np.max(y), vary = vary)
    #params.add('freq_offset', value = freq_offset, min = np.min(x), max = np.max(x), vary = vary)
    
    params.add('a', value = 1.0, vary = vary)
    params.add('b', value = 1.0, vary = vary)
    params.add('c', value = 1.0, vary = vary)
    
    
    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(x, y))
    result = minner.minimize()
    
    # Store the Confidence data from the fit
    con_report = lmfit.fit_report(result.params)
    
    (x_plot, model) = fcn2min(result.params, x, y, plot_fit = True)
    
    # get residuals
    (residuals) = fcn2min(result.params, x, y)
    
    #:print(result.params)
    
    return (x_plot, model, result, residuals)


# define objective function: returns the array to be minimized
# function to minimize
def fcn2min_single_exp(params, x, data, plot_fit = False):
    
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value

    if plot_fit == False:
        x_fit = np.array(x)
    else:
        no_points = 1000
        x_plot = np.linspace(np.min(x), np.max(x), no_points)
        x_fit = x_plot

    model = a * np.exp(-b * x_fit) + c

    if plot_fit == False:
        return model - data
    else:
        return (x_plot, model)



def fit_theory_single_exp(x, y, vary = True):
    
    # fits a sum of gaussians to a data set
    # my_lines is a list of frequency offsets

    params = Parameters()
    
    #params.add('a', value = np.min(y), min = 0*np.min(y), max = np.max(y), vary = vary)
    #params.add('freq_offset', value = freq_offset, min = np.min(x), max = np.max(x), vary = vary)
    
    if y[0] > y[-1]:
        a =  np.max(y)
        b =  1.0
        c =  0.0
    else:
        a = -np.max(y)
        b =  1.0
        c = np.max(y)

    params.add('a', value = a, vary = vary)
    params.add('b', value = b, vary = vary)
    params.add('c', value = c, vary = vary)
    
    
    # do fit, here with leastsq model
    minner = Minimizer(fcn2min_single_exp, params, fcn_args=(x, y))
    result = minner.minimize()
    
    # Store the Confidence data from the fit
    con_report = lmfit.fit_report(result.params)
    
    (x_plot, model) = fcn2min_single_exp(result.params, x, y, plot_fit = True)
    
    # get residuals
    (residuals) = fcn2min_single_exp(result.params, x, y)
    
    #:print(result.params)
    
    return (x_plot, model, result, residuals)


#################################################################################




data = [
        {'date':20210108, 'time':111054, 'title' : 'AlCl Scan'},
        {'date':20210122, 'time':131636, 'title' : 'AlCl Scan'},
        {'date':20210128, 'time':163613, 'title' : 'AlCl Scan'},
        {'date':20210203, 'time':151628, 'title' : 'AlCl Scan'},
        {'date':20210203, 'time':164112, 'title' : 'Al Scan'}, # Al scan
        {'date':20210210, 'time':125419, 'title' : 'AlCl Scan'}, # AlCl scan
        {'date':20210210, 'time':152249, 'title' : 'Al Scan'}, # Al scan
        {'date':20210210, 'time':164806, 'title' : 'K Scan'}, # K scan
        {'date':20210212, 'time':111519, 'title' : 'AlCl Scan'}, # AlCl scan
        {'date':20210212, 'time':125709, 'title' : 'Al Scan'}, # Al scan
        {'date':20210212, 'time':142438, 'title' : 'K Scan'}, # K scan
        
        {'date':20210218, 'time':103239, 'title' : 'AlCl Scan'}, # AlCl scan
        {'date':20210218, 'time':123517, 'title' : 'Al Scan'}, # Al scan
        {'date':20210218, 'time':140601, 'title' : 'K Scan'}, # K scan
    ]

options = { 'offset_avg_points' : 20, 'moving_avg_no' : 0, 'subtract_dc_offset' : True }


size_x = 0.075
size_y = 0.075

photodiode_offset_AlCl = 3.86e-3
photodiode_offset_Al_K = -6.42e-4

photodiode_offset_Al_K2 = -9.94e-4

photodiode_offset_Al_K3 = -1.38e-3

#t_int = [10.885, 11.00]
t_int = [10.9, 11.5]

result = []
for k in range(len(data)):
    result.append(get_img_data(data[k], options))


plot_options = [
        { 
            'photodiode_offset' : photodiode_offset_AlCl, 't_integrate' : t_int,  'channel' : 0, 
            'img_integrate' : [ 
                [ [5.6, 4.65], [size_x, size_y], 'Al/KCl, 1:4' ],
                [ [5.85, 5.1], [size_x, size_y], 'Al/MgCl2, 2:1' ],
                [ [5.4, 5.1], [size_x, size_y], 'AlCl3' ],
                [ [5.5, 5.4], [size_x, size_y], 'BG' ]
                ]
            },
        { 
            'photodiode_offset' : photodiode_offset_AlCl, 't_integrate' : t_int, 'channel' : 0, 
            'img_integrate' : [ 
                [ [5.864, 4.842], [size_x, size_y], 'Al/KCl, 1:3' ],
                [ [5.51, 5.067], [size_x, size_y], 'Al/KCl, 1:1' ],
                [ [5.478, 4.6], [size_x, size_y], 'Al/KCl, 3:1' ],
                [ [5.7, 5.3], [size_x, size_y], 'BG' ]
                ]
        },
        { 'photodiode_offset' : photodiode_offset_AlCl, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [5.8, 4.65], [size_x, size_y], 'Al/KCl, 1:10' ],
                        [ [5.4, 5.0], [size_x, size_y], 'Al/KCl, 1:3' ],
                        [ [5.85, 5.1], [size_x, size_y], 'Al/KCl, 0' ]
                        ]
        },
{ 'photodiode_offset' : photodiode_offset_AlCl, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [5.7, 4.7], [size_x, size_y], 'Al/KCl, 1:20' ],
                        [ [5.55, 4.3], [size_x, size_y], 'Al/KCl, 1:3' ],
                        [ [5.25, 4.65], [size_x, size_y], 'Al/KCl, 10:1' ],
                        [ [5.2, 4.25], [size_x, size_y], 'Al/KCl, 1:2' ]
                        ]
        },
{ 'photodiode_offset' : photodiode_offset_Al_K, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [5.7, 4.7], [size_x, size_y], 'Al/KCl, 1:20' ],
                        [ [5.55, 4.3], [size_x, size_y], 'Al/KCl, 1:3' ],
                        [ [5.25, 4.65], [size_x, size_y], 'Al/KCl, 10:1' ],
                        [ [5.2, 4.25], [size_x, size_y], 'Al/KCl, 1:2' ]
                        ]
        },

{ 'photodiode_offset' : photodiode_offset_AlCl, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [4.55, 4.62], [size_x, size_y], 'Al/KCl, 8:1' ],
                        [ [4.474, 4.2], [size_x, size_y], 'Al/KCl, 10:1' ],
                        [ [4.964, 4.163], [size_x, size_y], 'Al/KCl, 5:1' ],
                        [ [4.997, 4.669], [size_x, size_y], 'Al/KCl, 1:4' ]
                        ]
        }, 
{ 'photodiode_offset' : photodiode_offset_Al_K, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [4.55, 4.62], [size_x, size_y], 'Al/KCl, 8:1' ],
                        [ [4.474, 4.2], [size_x, size_y], 'Al/KCl, 10:1' ],
                        [ [4.964, 4.163], [size_x, size_y], 'Al/KCl, 5:1' ],
                        [ [4.997, 4.669], [size_x, size_y], 'Al/KCl, 1:4' ]
                        ]
        }, 
{ 'photodiode_offset' : photodiode_offset_Al_K, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [4.55, 4.62], [size_x, size_y], 'Al/KCl, 8:1' ],
                        [ [4.474, 4.2], [size_x, size_y], 'Al/KCl, 10:1' ],
                        [ [4.964, 4.163], [size_x, size_y], 'Al/KCl, 5:1' ],
                        [ [4.997, 4.669], [size_x, size_y], 'Al/KCl, 1:4' ]
                        ]
        }, 
{ 'photodiode_offset' : photodiode_offset_AlCl, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [4.55, 4.62], [size_x, size_y], 'Al/KCl, 1:5' ],
                        [ [4.55, 4.2], [size_x, size_y], 'Al/KCl, 1:8' ],
                        [ [4.9, 4.24], [size_x, size_y], 'KCl' ],
                        [ [4.997, 4.6], [size_x, size_y], 'Al/KCl, 1:3' ]
                        ]
        }, 
{ 'photodiode_offset' : photodiode_offset_Al_K2, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [4.55, 4.62], [size_x, size_y], 'Al/KCl, 1:5' ],
                        [ [4.55, 4.2], [size_x, size_y], 'Al/KCl, 1:8' ],
                        [ [4.9, 4.24], [size_x, size_y], 'KCl' ],
                        [ [4.997, 4.6], [size_x, size_y], 'Al/KCl, 1:3' ]
                        ]
        }, 
{ 'photodiode_offset' : photodiode_offset_Al_K2, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [4.55, 4.62], [size_x, size_y], 'Al/KCl, 1:5' ],
                        [ [4.55, 4.2], [size_x, size_y], 'Al/KCl, 1:8' ],
                        [ [4.9, 4.24], [size_x, size_y], 'KCl' ],
                        [ [4.997, 4.6], [size_x, size_y], 'Al/KCl, 1:3' ]
                        ]
        },

# 20210218
 { 'photodiode_offset' : photodiode_offset_AlCl, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [4.5, 4.55], [size_x, size_y], 'Al/KCl, 1:50' ],
                        [ [4.45, 4.2], [size_x, size_y], 'Al/KCl, 1:1' ],
                        [ [5.00, 4.6], [size_x, size_y], 'Al/KCl, 2:1' ],
                        [ [4.95, 4.24], [size_x, size_y], 'Al/KCl, 7:1' ]
                        ]
        }, 
{ 'photodiode_offset' : photodiode_offset_Al_K3, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [4.5, 4.55], [size_x, size_y], 'Al/KCl, 1:50' ],
                        [ [4.45, 4.2], [size_x, size_y], 'Al/KCl, 1:1' ],
                        [ [5.00, 4.6], [size_x, size_y], 'Al/KCl, 2:1' ],
                        [ [4.95, 4.24], [size_x, size_y], 'Al/KCl, 7:1' ]
                        ]
        },  
{ 'photodiode_offset' : photodiode_offset_Al_K3, 't_integrate' : t_int, 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [4.5, 4.55], [size_x, size_y], 'Al/KCl, 1:50' ],
                        [ [4.45, 4.2], [size_x, size_y], 'Al/KCl, 1:1' ],
                        [ [5.00, 4.6], [size_x, size_y], 'Al/KCl, 2:1' ],
                        [ [4.95, 4.24], [size_x, size_y], 'Al/KCl, 7:1' ]
                        ]
        },  
        ]

# [scan number, box number, Al part, KCl part]

AlCl_scans = [ 
        [0, 0, 1, 4],
        [1, 0, 1, 3],
        [1, 1, 1, 1],
        [1, 2, 3, 1],
        [2, 0, 1, 10],
        [2, 1, 1, 3],
        [2, 2, 0, 1], 
        [3, 0, 1, 20], 
        [3, 1, 1, 3],
        [3, 2, 10, 1], 
        [3, 3, 1, 2],
        [5, 0, 8, 1], 
        [5, 1, 10, 1],
        [5, 2, 5, 1], 
        [5, 3, 1, 4],
        [8, 0, 1, 5], 
        [8, 1, 1, 8],
        [8, 2, 0, 1], 
        [8, 3, 1, 3],
        [11, 0, 1, 50], 
        [11, 1, 1, 1],
        [11, 2, 2, 1], 
        [11, 3, 7, 1],
        ]

Al_scans = [ 
        [6, 0, 8, 1], 
        [6, 1, 10, 1],
        [6, 2, 5, 1], 
        [6, 3, 1, 4],
        [9, 0, 1, 5], 
        [9, 1, 1, 8],
        [9, 2, 0, 1], 
        [9, 3, 1, 3],
        [12, 0, 1, 50], 
        [12, 1, 1, 1],
        [12, 2, 2, 1], 
        [12, 3, 7, 1],
        ]

K_scans = [ 
        [7, 0, 8, 1], 
        [7, 1, 10, 1],
        [7, 2, 5, 1], 
        [7, 3, 1, 4],
        [10, 0, 1, 5], 
        [10, 1, 1, 8],
        [10, 2, 0, 1], 
        [10, 3, 1, 3],
        [13, 0, 1, 50], 
        [13, 1, 1, 1],
        [13, 2, 2, 1], 
        [13, 3, 7, 1],
        ]


yields = []

for k in range(len(result)):
    (img, abs_yields) = plot_target_img(result[k], plot_options[k], data[k]['title'], do_plot = True)

    yields.append(abs_yields)



x_AlCl = []
y_AlCl = []

for k in range(len(AlCl_scans)):

    scan_id = AlCl_scans[k][0]
    box_id = AlCl_scans[k][1]
    al_part = AlCl_scans[k][2]
    kcl_part = AlCl_scans[k][3]

    x_AlCl.append(np.float(al_part) / np.float(kcl_part))
    y_AlCl.append(yields[scan_id][box_id])

x_Al = []
y_Al = []

for k in range(len(Al_scans)):

    scan_id = Al_scans[k][0]
    box_id = Al_scans[k][1]
    al_part = Al_scans[k][2]
    kcl_part = Al_scans[k][3]

    x_Al.append(np.float(al_part) / np.float(kcl_part))
    y_Al.append(yields[scan_id][box_id])



x_K = []
y_K = []

for k in range(len(K_scans)):

    scan_id = K_scans[k][0]
    box_id = K_scans[k][1]
    al_part = K_scans[k][2]
    kcl_part = K_scans[k][3]

    x_K.append(np.float(al_part) / np.float(kcl_part))
    y_K.append(yields[scan_id][box_id])





# plot all yields

plt.figure(figsize = (10, 6))



plt.subplot(2,2,2)

plt.plot(x_Al, y_Al, 'bx')

(xfit, yfit, result, residuals) = fit_theory_single_exp(x_Al, y_Al)

Al_par = result.params

plt.plot(xfit, yfit, 'r')


plt.ylabel('Absorption yield Al (%)')
plt.xlabel('Al/Cl Molar Ratio ')

plt.ylim(0,110)

plt.tight_layout()


plt.subplot(2,2,3)
plt.plot(x_K, y_K, 'bx')

(xfit, yfit, result, residuals) = fit_theory_single_exp(x_K, y_K)

K_par = result.params

plt.plot(xfit, yfit, 'r')


plt.ylabel('Absorption yield K (%)')
plt.xlabel('Al/Cl Molar Ratio ')

plt.ylim(0,110)

plt.tight_layout()



plt.subplot(2,2,1)
plt.plot(x_AlCl[0:11], y_AlCl[0:11], 'bo')
plt.plot(x_AlCl[11:], y_AlCl[11:], 'rx')

(xfit, yfit, result, residuals) = fit_theory(x_AlCl, y_AlCl)
plt.plot(xfit, yfit, 'r')


#x = np.linspace(0, 10, 100)
#
#y  = K_par['a'].value * np.exp(-K_par['b'].value * x) + K_par['c'].value
#y *= Al_par['a'].value * np.exp(-Al_par['b'].value * x) + Al_par['c'].value
#
#plt.plot(x, y, 'k--')

plt.ylabel('Absorption yield AlCl (%)')
plt.xlabel('Al/Cl Molar Ratio ')

plt.ylim(0,110)

plt.tight_layout()



###########################################

plt.figure(figsize = (10, 6))

plt.subplot(2,1,1)

#print(AlCl_scans[11:15])

#plt.plot(y_Al, y_AlCl[11:15], 'bo', markersize = 10)

plt.scatter(y_Al[0:4],  y_AlCl[11:15], color = 'b', marker = 'o', fc = 'w', lw = 2.0, s = 40)
plt.scatter(y_Al[4:8],  y_AlCl[15:19], color = 'b', marker = 'o', fc = 'w', lw = 2.0, s = 40)
plt.scatter(y_Al[8:12], y_AlCl[19:23], color = 'b', marker = 'o', fc = 'w', lw = 2.0, s = 40)

plt.plot([0, 100], [100, 0], 'r--')

plt.xlabel('Absorption yield Al (%)')
plt.ylabel('Absorption yield AlCl (%)')

plt.xlim(0,110)
plt.ylim(0,110)

plt.tight_layout()

plt.subplot(2,1,2)

#plt.plot(y_K, y_AlCl[11:15], 'bo', markersize = 10)

plt.scatter(y_K[0:4], y_AlCl[11:15], color = 'b', marker = 'o', fc = 'w', lw = 2.0, s = 40)
plt.scatter(y_K[4:8], y_AlCl[15:19], color = 'b', marker = 'o', fc = 'w', lw = 2.0, s = 40)
plt.scatter(y_K[8:12], y_AlCl[19:23], color = 'b', marker = 'o', fc = 'w', lw = 2.0, s = 40)


plt.plot([0, 100], [0, 100], 'r--')

plt.xlabel('Absorption yield K (%)')
plt.ylabel('Absorption yield AlCl (%)')

plt.xlim(0,110)
plt.ylim(0,110)

plt.tight_layout()




plt.show()



