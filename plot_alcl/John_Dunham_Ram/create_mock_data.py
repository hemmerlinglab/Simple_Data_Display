import numpy as np
import matplotlib.pyplot as plt


def get_E(Y,v,J,mat_size):
    #print(Y)
    E = 0.0
    #d = 4
    for i in range(mat_size[0]):
        for j in range(mat_size[1]):
            E += Y[i][j] * (v+0.5)**i * (J*(J+1.0))**j

    return E


# ground state
varr = [0,1]
Jarr = np.arange(0, 15)

# excited state
v2arr = [0,1]
J2arr = np.arange(0, 15)

Yg00 = 0.0
Yg10 = 481.774655196
Yg01 = 2.4393006612e-1
Yg11 = -1.6110822121e-3

Ye00 = 38237.0
Ye10 = 570.774655196
Ye01 = 2.3393006612e-1
Ye11 = -1.6110822121e-3



Yg = np.array([[Yg00, Yg01], [Yg10, Yg11]])
Ye = np.array([[Ye00, Ye01], [Ye10, Ye11]])


mat_size = [2,2]

dE_arr = []
Ps = []
Qs = []
Rs = []
qs = [[],[],[],[]]

for v in varr:
    for v2 in v2arr:
        for J in Jarr:
            for J2 in J2arr:

                if np.abs(J2 - J) <= 1:

                    dE = get_E(Ye, v2, J2, mat_size) - get_E(Yg, v, J, mat_size)

                    dE_arr.append(dE)

                    qs[0].append(v2)
                    qs[1].append(J2)
                    qs[2].append(v)
                    qs[3].append(J)

                    if J2 == J:
                        Qs.append([J, dE])
                    if J2 - 1 == J:
                        Ps.append([J, dE])
                    if J2 + 1 == J:
                        Rs.append([J, dE])
                        
Ps = np.array(Ps)
Qs = np.array(Qs)
Rs = np.array(Rs)


# plot lines
plt.figure()
plt.plot(Ps[:, 0], Ps[:, 1], 'bx')
plt.plot(Qs[:, 0], Qs[:, 1], 'rx')
plt.plot(Rs[:, 0], Rs[:, 1], 'kx')


# test script

mat_size = [2, 2]

from ram_dunham_fit import *

result_params, ps, cr = fit_dunham(qs,dE_arr, mat_size = mat_size)

#print(rps)
#print(ps)
print(cr)


fit_y = fcn2min(result_params, qs, dE_arr, mat_size = mat_size, get_fit = True)

print(fit_y)

fit_x = qs[1] # all gound state J's

plt.plot(fit_x, fit_y, 'o')

plt.show()



