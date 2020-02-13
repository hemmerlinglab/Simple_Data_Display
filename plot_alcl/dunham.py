from fit_dunham import *
import numpy as np

c = 299792458


def print_params(par, state = 'Yg'):

    print("Te   = {0:12.6f} THz".format(100*c*par[state + '00'].value/1e12))
    print("we   = {0:12.6f} THz".format(100*c*par[state + '10'].value/1e12))
    print()
    print("wexe = {0:12.4f} GHz".format(100*c*par[state + '20'].value/1e9))
    print("weye = {0:12.4f} GHz".format(100*c*par[state + '30'].value/1e9))
    print("weze = {0:12.4f} GHz".format(100*c*par[state + '40'].value/1e9))
    
    print()
    print("Be   = {0:12.4f} GHz".format(100*c*par[state + '01'].value/1e9))
    print("ae   = {0:12.4f} GHz".format(100*c*par[state + '11'].value/1e9))




# d = [ [gs_v, gs_J, ex_v, ex_J, freq] ]

d = [ 
        [0,0,0,0, 38237.00],
        [1,0,1,0, 38201.72],
        [2,0,2,0, 38157.21],
        [3,0,3,0, 38105.26],
        [4,0,4,0, 38043.44],
        [4,0,5,0, 38433.85],
        [5,0,5,0, 37972.12], # typo
        [6,0,6,0, 37888.98],
        [7,0,6,0, 37435.88], 
        [9,0,7,0, 36898.34],
        [10,0,8,0, 36794.91],
        [15,0,10,0, 35248.74]
        ]

# rotational states

d.append([4,40,4,40,38036.01])
d.append([4,48,4,48,38032.55])
d.append([6,30,6,30,37880.80])
d.append([9,25,7,25,36892.90])




(result, guesses) = do_fit(d)

print()
for k in result.keys():

    print("{0} = {1:6.4f}".format(k, result[k].value))

print()

# print results

total_err = 0

for n in range(len(d)):

    fit = fcn2min(result, [], [d[n]], get_fit = True)[0]

    err = d[n][4] - fit

    total_err += np.abs(err)

    print('(' + "({0},{1}) -> ({2},{3})".format(d[n][0], d[n][1], d[n][2], d[n][3]) + ') = ' + "{0:8.2f} / {1:8.2f}   diff: {2:8.2f}".format(d[n][4], fit, err))


print("\nTotal error = {0:6.2f}\n".format(total_err))




# values of fit in MHz

print("*"*40)
print("Fit results")
print("*"*40)
print()

print("*"*40)
print("Ground state")
print("*"*40)
print_params(result, 'Yg')

print()
print("*"*40)
print("Excited state")
print("*"*40)
print_params(result, 'Ye')




print()


