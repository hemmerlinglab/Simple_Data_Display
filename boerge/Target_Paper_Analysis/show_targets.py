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
        #{'date':20210108, 'time':111054},
        {'date':20210122, 'time':122819},
        {'date':20210122, 'time':131636}
    ]

options = { 'offset_avg_points' : 20, 'moving_avg_no' : 0, 'subtract_dc_offset' : True }


size_x = 0.05
size_y = 0.05

photodiode_offset = 3.86e-3

result = []
for k in range(len(data)):
    result.append(get_img_data(data[k], options))


plot_options = { 'photodiode_offset' : photodiode_offset, 't_integrate' : [10.885, 11.00], 'channel' : 0, 
        'img_integrate' : [ 
            #[   [5.7, 4.65], [0.125, 0.05]  ],
            #[   [5.85, 5.1], [0.125, 0.05]  ],
            #[   [5.45, 5.2], [0.05, 0.125]  ]
            [ [5.864, 4.842], [size_x, size_y] ],
            [ [5.51, 5.067], [size_x, size_y] ],
            [ [5.478, 4.6], [size_x, size_y] ],
            [ [5.7, 5.3], [size_x, size_y] ],
        ]}

for k in range(len(result)):
    plot_target_img(result[k], plot_options)


plt.show()



