import numpy as np
from aux_functions import get_dunham, get_energy


c = 299792458
kB = 1.3806482e-23
h_planck = 6.63607e-34
amu = 1.66e-27

massAl = 26.98153841
massCl_35 = 34.96885269
massCl_37 = 36.96590258


def get_pop(v, J, we, B, T):

    # returns population distribution

    Ewe = we * h_planck
    EB = B * h_planck

    beta = 1/(kB*T)

    Ev = Ewe * (v+0.5)
    
    Erot = lambda J : EB * J * (J+1)

    # rotational distribution
    J_arr = np.arange(0, 200, 1)

    #print(Jarr)

    Nrot = np.sum( (2*J_arr+1) * np.exp(-beta * Erot(J_arr)) )
    
    pop_rot = 1/Nrot * (2*J+1) * np.exp(-beta * Erot(J))

    # vibrational distribution
    Nvib = np.exp(beta * Ewe/2)/(np.exp(beta * Ewe) - 1)

    pop_vib = 1/Nvib * np.exp(-beta * Ev)

    return pop_vib * pop_rot



def get_doppler(T, f0):
    # returns Doppler width for 27Al35Cl

    return np.sqrt(8*kB*T*np.log(2)/((35+27)*amu*c**2)) * f0



def simple_gauss(x, x0, A, w):
    return A * np.exp( -(x-x0)**2/(2*w**2) )


def get_gaussians(x, cnts, a, T = 4):

    y = np.zeros(len(x))

    for k in range(len(cnts)):

         w = get_doppler(T, cnts[k])
         y += simple_gauss(x, cnts[k], a[k], w)

    return y



def get_line_centers(Yg, Ye, ve, vg, Jmax = 1, T = 50, real_amplitudes = True):

    Jg_arr = np.arange(0, Jmax+1)    
    Je_arr = np.arange(1, Jmax+1)
    
    f_P = []
    f_Q = []
    f_R = []
    
    A_P = []
    A_Q = []
    A_R = []
    
    for Jg in Jg_arr:
        for Je in Je_arr:
            
            # only dipole transitions
            eng = get_energy(Yg, Ye, vg, Jg, ve, Je)
            
            # apply population of ground states
            if real_amplitudes:
                A = get_pop(vg, Jg, 100*c*Yg[1][0], 100*c*Yg[0][1], T)
            else:
                A = 1.0

            if Je - Jg == -1: 
                f_P.append(eng)
                A_P.append(A)
    
            if Je - Jg ==  0: 
                f_Q.append(eng)
                A_Q.append(A)
    
            if Je - Jg == +1: 
                f_R.append(eng)
                A_R.append(A)
    
    f_P = np.array(f_P)
    f_Q = np.array(f_Q)
    f_R = np.array(f_R)
    
    A_P = np.array(A_P)
    A_Q = np.array(A_Q)
    A_R = np.array(A_R)

    return ([f_P, f_Q, f_R], [A_P, A_Q, A_R])



def get_spectrum(Ug = [], Ue = [], ve = 0, vg = 0, Jmax = 10, T = 10, df = 100e9, no_of_points = 5000):

    (Yg35, Ye35, Yg37, Ye37) = get_dunham(Ug, Ue)
   
    (lines35, ampl35) = get_line_centers(Yg35, Ye35, ve, vg, Jmax = Jmax, T = T)
    (lines37, ampl37) = get_line_centers(Yg37, Ye37, ve, vg, Jmax = Jmax, T = T)

    cnt_freq = np.mean(lines35[1])

    nus = np.linspace(cnt_freq - df, cnt_freq + df, no_of_points)
    
    g35 = []
    g37 = []
    for k in [0,1,2]:
        g35.append(get_gaussians(nus, lines35[k], 0.76*ampl35[k], T = T/2))
        g37.append(get_gaussians(nus, lines37[k], 0.24*ampl37[k], T = T/2))

    return (nus, g35, g37)


