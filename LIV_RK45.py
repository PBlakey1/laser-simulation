import numpy as np
import matplotlib.pyplot as plt
## Define Variables
xi = 0.2                #Field Confinment Factor
V_a = 7.5e-17             #Volume of Active region
e = 1.6e-19             #Charge of an Electron in Coulombs
c = 2.998e8             #Speed of light in m/s
n_eq = 3.59              #Effective index of lasing mode.
a = 2.75e-12            #Tangential Coefficent
L = 0.0003               #Cavity length (meters)
Rb = 0.9                #First Mirror reflectivity  (Intensity)
Rf = 0.9                #Second Mirror reflectivity (Intensity)
alpha_loss = 1/8.68     #Round trip loss coefficent
nu = c/n_eq             #Propagation Speed of optical wave in laser (speed of light / effective index)
N_g = 1.89e8/V_a            #Transparent electron density
alpha = 2
T = 298.2
N_c = 3.65e18
tau_s = 2.79e-9

hbar = 1.055e-34
omega = 2*np.pi*3.527e14

t_span = (0,20e-9)

dt = tau_s/20000000



#G_th_1 = nu*(alpha_loss +1/(2*L)*np.log(1/(Rb*Rf)))
G_th = 5.01e11
N_th = G_th/(xi*a) + N_g
I_th = e*N_th*V_a/tau_s
I = 1.5*I_th
## Simulation
## The following simulation attempts to solve the following rate eqations for Photon number and electron density
## dS_p/dt = (G_p - G_th)S_p + xi*a*N. Note xi*a*N is an approximation of the spontaneous emission term
## dN/dt = -G_p S_p/V_a - N/tau_s + I/(e*V_a)
## The gain G_p = xi*a((N-N_g)

def rk45_step(fun, t, y, h, *args):
    k1 = h * fun(t, y, *args)
    k2 = h * fun(t + 0.25 * h, y + 0.25 * k1, *args)
    k3 = h * fun(t + (3 / 8) * h, y + (3 / 32) * k1 + (9 / 32) * k2, *args)
    k4 = h * fun(t + 12/13 * h, y + (1932/2197) * k1 - 7200/2197 * k2 + 7296/2197 * k3, *args)
    k5 = h * fun(t + h, y + 439/216 * k1 - 8 * k2 + 3680/513 * k3 - 845/4104 * k4, *args)
    k6 = h * fun(t + 0.5 * h, y - 8/27 * k1 + 2 * k2 + (3544/2565) * k3 +
                (1859/4104) * k4 - (11/40) * k5, *args)

    y_fourth = y + (25/216) * k1 + (1408/2565) * k3 + (2197/4101) * k4 - (1/5) * k5
    y_fifth = y + (16/135) * k1 + (6656/12825) * k3 + (28561/56430) * k4 - (9/50) * k5 + (2/55) * k6
    #error = np.sqrt((xi*a*(y_fifth[0] - y_fourth[0]))**2 + (y_fifth[1] - y_fourth[1])**2)
    error = np.sqrt((2*(y_fifth[0] - y_fourth[0])/(y_fifth[0] + y_fourth[0])) ** 2 + (2*(y_fifth[1] - y_fourth[1])/(y_fifth[1] + y_fourth[1])) ** 2)
    return y_fourth, error

def rate_equations(t,y,xi,a,N_g,G_th,V_a,tau_s,I,e):
    N, S = y
    #term_1 = -xi*a*(N - N_g) * S / V_as
    #term_2 = - N / tau_s
    #term_3 =I / (e * V_a)
    #dN_dt =  -xi*a*(N - N_g) * S / V_a - N / tau_s + I / (e * V_a)
    #dS_dt = (xi*a*(N - N_g) - G_th) * S + xi * a * N
    # dN_dt =  -xi*a*(N - N_g) * S / V_a - N / tau_s + I / (e * V_a)

    dN_dt2 = -xi*a*(N/V_a - N_g)*S - N/tau_s + I/e
    dS_dt = (xi * a * (N/V_a - N_g) - G_th) * S + xi * a * N/V_a
    return np.array([dN_dt2, dS_dt])
def solve_rk45(fun, t_span, y0, h, *args):
    t_values = [t_span[0]]
    y_values = [y0]
    y_values[0] = y0
    tol_h = 1e-6
    while t_values[-1] < t_span[1]:
        t = t_values[-1]
        y = y_values[-1]
        y_next, error = rk45_step(fun, t, y, h, *args)
        if error == 0:
            s = 1
        else:
            s = 0.84*(tol_h/error)**0.25
        t_values.append(t + h)
        y_values.append(y_next)
        h = s*h
    return np.array(t_values), np.array(y_values).T

# Initial conditions
#N0 = 0.0  # Initial carrier density
#S0 = 0.0  # Initial photon number
# Initial conditions
N0 = N_th*V_a # Initial carrier density
#S0 = 8.8e2  # Initial photon number
S0 = 0  # Initial photon number
# Solve the rate equations using RK45 method
t_values, y_values = solve_rk45(rate_equations, t_span, [N0, S0],dt,xi,a,N_g,G_th,V_a,tau_s,I,e)

N = y_values[0]
S_p = y_values[1]

plt.subplot(2, 1, 1)
plt.plot(t_values,S_p)
plt.subplot(2, 1, 2)
plt.plot(t_values, N)
plt.plot(t_values,V_a*N_th*np.ones(np.size(N)))
plt.show()
