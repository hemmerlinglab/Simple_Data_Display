

from fit_data import *
import numpy as np


d = [ 
        [0,0, 38237.00],
        [1,1, 38201.72],
        [2,2, 38157.21],
        [3,3, 38105.26],
        [4,4, 38043.44],
        [5,4, 38433.85],
        #[5,5, 38972.12],
        [6,6, 37888.98],
        [6,7, 37435.88],
        [7,9, 36898.34],
        [8,10, 36794.91],
        [10,15, 35248.74]
        ]

x = list(map(lambda s : [s[0], s[1]], d))
y = list(map(lambda s : s[2], d))


model = do_fit(x, y)


result = model.params


for k in result.keys():

    print(k + ' = ' + str(result[k].value))


# print results

for n in range(len(d)):

    fit = fcn2min(result, [x[n]], y[n], get_fit = True)[0]

    print('(' + str(d[n][0]) + ' -> ' + str(d[n][1]) + ') = ' + "{0:8.2f} / {1:8.2f}   diff: {2:8.2f}".format(d[n][2], fit, d[n][2] - fit))


