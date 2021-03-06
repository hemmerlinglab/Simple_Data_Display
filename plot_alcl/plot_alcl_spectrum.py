import numpy as np
import matplotlib.pyplot as plt
from fit_dunham import energy

c = 299792458
kB = 1.3806482e-23
h_planck = 6.63607e-34
amu = 1.66e-27

def get_pop(J, B, T, N0 = 1):

    pop = N0 * (2*J+1) * np.exp(-B * h_planck * J*(J+1)/(kB*T))

    return pop / (np.sqrt(kB*T/(2*h_planck*B)) - 1.0/2.0)

def get_doppler(T, f0):

    return np.sqrt(8*kB*T*np.log(2)/((35+27)*amu*c**2)) * f0


def gauss(x, x0, A, w):
    return A * np.exp( -(x-x0)**2/(2*w**2) )


def get_transitions(Yg, Ye, ve, vg, cnt_freqs, df = 100e6, T = 1, Jmax = 1):

    Jg_arr = np.arange(0, Jmax+1)    
    Je_arr = np.arange(0, Jmax+1)
    
    w = get_doppler(T, 100*c*cnt_freq)
    
    nus = np.linspace(100*c*(cnt_freq) - df, 100*c*(cnt_freq) + df, 2000)
    
    spectrum_P = np.zeros(len(nus))
    spectrum_Q = np.zeros(len(nus))
    spectrum_R = np.zeros(len(nus))
    
    f_P = []
    f_Q = []
    f_R = []
    
    for Jg in Jg_arr:
        for Je in Je_arr:
            
            # only dipole transitions
    
            eng = 100*c*(energy(Ye, ve, Je) - energy(Yg, vg, Jg))
    
            # apply population of ground states
            A = get_pop(Jg, 100*c*Yg[1][0], T)
    
            if Je - Jg == -1: 
                spectrum_P += gauss(nus, eng, A, w)
                f_P.append(eng)
                #plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'r' )
    
            if Je - Jg ==  0: 
                spectrum_Q += gauss(nus, eng, A, w)
                f_Q.append(eng)
                #plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'g' )
    
            if Je - Jg == +1: 
                spectrum_R += gauss(nus, eng, A, w)
                f_R.append(eng)
                #plt.plot( 2 * [(eng - 100*c*cnt_freq) / 1e9], [-0.1,0.0], 'b' )
    
    
    f_P = np.array(f_P) - 100*c*cnt_freq
    f_Q = np.array(f_Q) - 100*c*cnt_freq
    f_R = np.array(f_R) - 100*c*cnt_freq

    nus = nus - 100*c*cnt_freqs

    return (nus, [spectrum_P, spectrum_Q, spectrum_R], [f_P, f_Q, f_R])


def plot_spectrum(nus, s, ve = 0, vg = 0, T = 0, cnt_freq = 0, style = '-', txt = '', abundance = 1.0):
    
    plt.plot( nus / 1e9, abundance * s[0], 'r' + style, label = 'P' + txt)
    plt.plot( nus / 1e9, abundance * s[1], 'g' + style, label = 'Q' + txt)
    plt.plot( nus / 1e9, abundance * s[2], 'b' + style, label = 'R' + txt)

    plt.xlabel("Frequency (GHz) - {0:6.6f} THz".format(100*c*cnt_freq/1e12))

    plt.title("AlCl transitions @ T = {2:2} K for v = {0} -> v' = {1} (J->J')".format(vg, ve, T))

    plt.xlim(np.min(nus)/1e9, np.max(nus)/1e9)

    plt.ylim(-0.1, 3.0)

    plt.legend()

def plot_transitions(f, cut = None, style = '-', txt = ''):
    
    if cut == None:
       cut = len(f[0])

    plt.plot(f[0][0:cut]/1e9, 'ro' + style, label = 'P' + txt)
    plt.plot(f[1][0:cut]/1e9, 'gx' + style, label = 'Q' + txt)
    plt.plot(f[2][0:cut]/1e9, 'bd' + style, label = 'R' + txt)

    plt.xlabel('Rotational number J')
    plt.ylabel("Frequency (GHz) - {0:6.6f} THz".format(100*c*cnt_freq/1e12))

    plt.legend(loc = 'upper right')



def scale_coeff(Y, mass1, mass2, scale = True):

    # scale coefficients with the reduced mass
    mu = (mass1 * mass2)/(mass1 + mass2)

    nvib = len(Y[0])

    if scale == False:

        # return U_kl = 1/mu^(-k/2 - l) * Y_kl
        U = []
        hlp = []
        for k in range(nvib):
            l = 0
            hlp.append( 1/(mu**(-k/2 - l)) * Y[0][k] )
        U.append(hlp)

        hlp = []
        for k in [0,1]:
            l = 1
            hlp.append( 1/(mu**(-k/2 - l)) * Y[1][k] )
        U.append(hlp)

        return U
    else:

        # return Y_kl = mu^(-k/2 - l) * U_kl

        U = Y.copy()
        Y = []
        hlp = []
        for k in range(nvib):
            l = 0
            hlp.append( (mu**(-k/2 - l)) * U[0][k] )
        Y.append(hlp)

        hlp = []
        for k in [0,1]:
            l = 1
            hlp.append( (mu**(-k/2 - l)) * U[1][k] )
        Y.append(hlp)

        return Y




###########################################################################################################################################################

# AlCl constants

#Yg = [[0.0, 482.0794801498224, -2.0285504669476797, 4.2244382524958823e-07, -1.4200906251997294e-05], [0.23086415576356475, -0.0015452837298277023]]
Ye = [[38250.600844308836, 457.2272074915964, -7.519552217461303, 0.29489258632225485, -0.02723237768418022], [0.23713820022739715, -0.003944286505039427]]

# from Bernath
Yg = [[0.0, 481.774655, -2.101, 6.638e-3, -2.025e-5],[0.2439300, -1.611e-3]]

massAl = 26.98153841
massCl_35 = 34.96885269
massCl_37 = 36.96590258

Ug = scale_coeff(Yg, massAl, massCl_35, scale = False)
Ue = scale_coeff(Ye, massAl, massCl_35, scale = False)



Yg35 = scale_coeff(Ug, massAl, massCl_35, scale = True)
Ye35 = scale_coeff(Ue, massAl, massCl_35, scale = True)

Yg37 = scale_coeff(Ug, massAl, massCl_37, scale = True)
Ye37 = scale_coeff(Ue, massAl, massCl_37, scale = True)






Jmax = 30
T = 10 # Kelvin
cnt_freq = 38237.0 # 1/cm
df = 100e9 # HZ



ve = 0
vg = 1
(nus10_35, spectrum10_35, f_lines_10_35) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, df, T, Jmax)
print((f_lines_10_35[1] + 100*c*cnt_freq)/3e12) # Q transition


plt.figure()
plot_spectrum(nus10_35, spectrum10_35, ve, vg, T, 0*cnt_freq, style = '-', txt = ' Cl35', abundance = 0.76)


# vibrational states
ve = 0
vg = 0
(nus35, spectrum35, f_lines35) = get_transitions(Yg35, Ye35, ve, vg, cnt_freq, df, T, Jmax)
(nus37, spectrum37, f_lines37) = get_transitions(Yg37, Ye37, ve, vg, cnt_freq, df, T, Jmax)


plt.figure()
plot_spectrum(nus35, spectrum35, ve, vg, T, cnt_freq, style = '-', txt = ' Cl35', abundance = 0.76)
plot_spectrum(nus37, spectrum37, ve, vg, T, cnt_freq, style = '--', txt = ' Cl37', abundance = 0.24)

plt.figure()
plot_transitions(f_lines35, cut = 12, style = '-', txt = ' Cl35')
plot_transitions(f_lines37, cut = 12, style = '--', txt = ' Cl37')



exp_freqs = np.array([
382.081350,
382.086184,
382.090879,
382.093449,
382.095838,
382.098180, # Q00
382.10050, 
382.11035, 
382.11240, 
382.11510, 
382.11710, 
382.12000, 
382.12185, 
382.12500, 
382.12665, 
382.12990, 
382.13145, 
382.13480,
382.13550,
382.13630,
382.13980,
382.14115, 
382.14480, 
382.10285 
])


exp_freqs = 3*exp_freqs*1e12
exp_freqs = (exp_freqs - 100*c*cnt_freq)
exp_freqs = exp_freqs


# search for closest lines
Jassign = []
theory_lines = []
theory_lines.extend(f_lines35[0])
theory_lines.extend(f_lines35[1])
theory_lines.extend(f_lines35[2])
theory_lines.extend(f_lines37[0])
theory_lines.extend(f_lines37[1])
theory_lines.extend(f_lines37[2])

theory_J = []
theory_J.extend(np.arange(0,Jmax))
theory_J.extend(np.arange(0,Jmax+1))
theory_J.extend(np.arange(0,Jmax))

theory_J.extend(np.arange(0,Jmax))
theory_J.extend(np.arange(0,Jmax+1))
theory_J.extend(np.arange(0,Jmax))

exp_J = []

for k in range(len(exp_freqs)):

    hlp = np.abs(exp_freqs[k] - theory_lines)

    ind = np.where( np.min(hlp) == hlp )[0][0]

    exp_J.append(theory_J[ind])


plt.plot(exp_J, exp_freqs/1e9, 'kx', markersize = 15)

plt.xlabel('Rotational number J')
plt.ylabel("Frequency (GHz) - {0:6.6f} THz".format(100*c*cnt_freq/1e12))











plt.show()




