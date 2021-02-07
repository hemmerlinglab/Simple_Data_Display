import sys
import matplotlib.pyplot as plt
import matplotlib

sys.path.append("../Analysis_Scripts/")
from data_functions import *

font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)



data = [
        {'date':20210108, 'time':111054},
        #{'date':20210122, 'time':122819},
        {'date':20210122, 'time':131636},
        {'date':20210128, 'time':163613},
        {'date':20210203, 'time':151628},
        {'date':20210203, 'time':164112},
    ]

options = { 'offset_avg_points' : 20, 'moving_avg_no' : 0, 'subtract_dc_offset' : True }


size_x = 0.05
size_y = 0.05

photodiode_offset = 3.86e-3

result = []
for k in range(len(data)):
    result.append(get_img_data(data[k], options))


plot_options = [
        { 
            'photodiode_offset' : photodiode_offset, 't_integrate' : [10.885, 11.00], 'channel' : 0, 
            'img_integrate' : [ 
                [ [5.6, 4.65], [size_x, size_y], 'Al/KCl, 1:4' ],
                [ [5.85, 5.1], [size_x, size_y], 'Al/MgCl2, 2:1' ],
                [ [5.4, 5.1], [size_x, size_y], 'AlCl3' ],
                [ [5.5, 5.4], [size_x, size_y], 'BG' ]
                ]
            },
        { 
            'photodiode_offset' : photodiode_offset, 't_integrate' : [10.885, 11.00], 'channel' : 0, 
            'img_integrate' : [ 
                [ [5.864, 4.842], [size_x, size_y], 'Al/KCl, 1:3' ],
                [ [5.51, 5.067], [size_x, size_y], 'Al/KCl, 1:1' ],
                [ [5.478, 4.6], [size_x, size_y], 'Al/KCl, 3:1' ],
                [ [5.7, 5.3], [size_x, size_y], 'BG' ]
                ]
        },
        { 'photodiode_offset' : photodiode_offset, 't_integrate' : [10.885, 11.00], 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [5.8, 4.65], [size_x, size_y], 'Al/KCl, 1:10' ],
                        [ [5.4, 5.0], [size_x, size_y], 'Al/KCl, 1:3' ],
                        [ [5.85, 5.1], [size_x, size_y], 'Al/KCl, 0' ]
                        ]
        },
{ 'photodiode_offset' : photodiode_offset, 't_integrate' : [10.885, 11.00], 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [5.7, 4.7], [size_x, size_y], 'Al/KCl, 1:20' ],
                        [ [5.55, 4.3], [size_x, size_y], 'Al/KCl, 1:3' ],
                        [ [5.25, 4.65], [size_x, size_y], 'Al/KCl, 10:1' ],
                        [ [5.2, 4.25], [size_x, size_y], 'Al/KCl, 1:2' ]
                        ]
        },
{ 'photodiode_offset' : photodiode_offset, 't_integrate' : [10.885, 11.00], 'channel' : 0, 
                    'img_integrate' : [ 
                        [ [5.7, 4.7], [size_x, size_y], 'Al/KCl, 1:20' ],
                        [ [5.55, 4.3], [size_x, size_y], 'Al/KCl, 1:3' ],
                        [ [5.25, 4.65], [size_x, size_y], 'Al/KCl, 10:1' ],
                        [ [5.2, 4.25], [size_x, size_y], 'Al/KCl, 1:2' ]
                        ]
        }, 
        ]

# [scan number, box number, Al part, KCl part]

KCl_scans = [ 
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
        [3, 3, 1, 2] 
        ]


yields = []

for k in range(len(result)):
    (img, abs_yields) = plot_target_img(result[k], plot_options[k])

    yields.append(abs_yields)


# plot all yields

plt.figure()


x = []
y = []

for k in range(len(KCl_scans)):

    scan_id = KCl_scans[k][0]
    box_id = KCl_scans[k][1]
    al_part = KCl_scans[k][2]
    kcl_part = KCl_scans[k][3]

    x.append(np.float(al_part) / np.float(kcl_part))
    y.append(yields[scan_id][box_id])


plt.plot(x, y, 'o')

plt.ylabel('Absorption yield (%)')
plt.xlabel('Al/Cl Molar Ratio ')

a = 20

#xf = np.linspace(0, 3, 100)
#yf = 100*(a*xf)/(1+a*xf)**2
#
#plt.plot(xf, yf)


plt.show()



