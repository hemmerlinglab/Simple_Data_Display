import sys
import matplotlib.pyplot as plt
import matplotlib

sys.path.append("../Analysis_Scripts/")
from data_functions import *

font = {#'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)



data = {
    'date':20210113, 'time':151505,
    }

options = { 'offset_avg_points' : 5, 'moving_avg_no' : 0, 'subtract_dc_offset' : False }

result = get_img_data(data, options)



ch0 = result['channels'][0]


plt.plot(ch0[0, :])

plt.show()


offset_avg = np.mean(ch0)
offset_dev = np.std(ch0)


print("Photodiode offset {0:.2e} +/- {1:.2e}".format(offset_avg, offset_dev))




