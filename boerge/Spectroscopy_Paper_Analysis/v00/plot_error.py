import numpy as np
import pickle
import matplotlib.pyplot as plt



# read in data file
arr = pickle.load( open( 'error.pickle', "rb" ) )

set_f = arr['set_f']
act_f= arr['act_f']
act_avg_f = arr['act_avg_f']


x = []
y = []

avg_err = []
std_err = []

for k in range(len(set_f)):

    x.extend(set_f[k])
    y.extend(act_avg_f[k])

    avg_err.extend( [np.mean( set_f[k] - act_avg_f[k] )] )
    std_err.extend(  [np.std( set_f[k] - act_avg_f[k] )] )

x = np.array(x)
y = np.array(y)

err = x - y


avg_err = np.array(avg_err)
std_err = np.array(std_err)

ind = np.where(np.abs(err) < 75e6)

x = x[ind]
err = err[ind]


total_std_err = np.std(err)

plt.figure()

plt.plot(err/1e6)

plt.xlabel('Set Frequency + {0:2.6f} (MHz)'.format(np.mean(x)/1e12))
plt.ylabel('Set - Act Frequency (MHz)')


plt.figure()

plt.subplot(2,1,1)
plt.plot(avg_err/1e6)
plt.subplot(2,1,2)
plt.plot(std_err/1e6)

plt.axhline(np.mean(std_err)/1e6, ls = '--', color = 'r')
plt.axhline(total_std_err/1e6, ls = '--', color = 'r')

plt.xlabel('Set Frequency + {0:2.6f} (MHz)'.format(np.mean(x)/1e12))
plt.ylabel('Std of errors (MHz)')

print("Std error {0:2.6f} MHz".format(total_std_err/1e6))


plt.show()




