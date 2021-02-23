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
        {'date':20210113, 'time':151505},
        {'date':20210128, 'time':163125},
        {'date':20210203, 'time':172431},
        {'date':20210212, 'time':141519},
        {'date':20210212, 'time':141708},
        {'date':20210218, 'time':114636},
       ]

options = { 'offset_avg_points' : 5, 'moving_avg_no' : 0, 'subtract_dc_offset' : False }

result = []
for k in range(len(data)):
    result.append(get_img_data(data[k], options))


for k in range(len(result)):
    ch0 = result[k]['channels'][0]


    plt.plot(ch0[0, :])



    offset_avg = np.mean(ch0)
    offset_dev = np.std(ch0)


    print("Scan {3}_{2} : Photodiode offset {0:.2e} +/- {1:.2e}".format(offset_avg, offset_dev, data[k]['time'], data[k]['date']))


plt.show()


