import numpy as np
import matplotlib.pyplot as plt
## Define Variables
xi = 0.2                #Field Confinment Factor
V_a = 1e-16             #Volume of Active region
e = 1.6e-19             #Charge of an Electron in Coulombs
c = 2.998e8             #Speed of light in m/s
n_eq = 3.1              #Effective index of lasing mode.
a = 2.75e-12            #Tangential Coefficent
L = 0.001               #Cavity length (meters)
Rb = 0.9               #First Mirror reflectivity  (Intensity)
Rf = 0.9               #Second Mirror reflectivity (Intensity)
alpha_loss = 1/8.68     # Round trip loss coefficent
nu = c/n_eq             #Propagation Speed of optical wave in laser (speed of light / effective index)
N_g = 2.1e24            #Transparent electron density
T = 298.2
N_c = 3.65e18;
tau_s = 1e3*2.85e7*T/(N_c*1e6)


hbar = 1.055e-34
omega = 2*np.pi*3.527e14

t_span = 1e-8
n_steps = 2**21

dt = t_span/n_steps
dt_next = dt
time = np.zeros(n_steps)
S_p = np.zeros(n_steps)
G_p = np.zeros(n_steps)
step_size = np.zeros(n_steps)
N = np.zeros(n_steps)

G_th = nu*(alpha_loss +1/(2*L)*np.log(1/(Rb*Rf)))
N_th = G_th/(xi*a) + N_g
I_th = e*N_th*V_a/tau_s
I = 1.5*I_th
## Simulation
## The following simulation attempts to solve the following rate eqations for Photon number and electron density
## dS_p/dt = (G_p - G_th)S_p + xi*a*N. Note xi*a*N is an approximation of the spontaneous emission term
## dN/dt = -G_p S_p/V_a - N/tau_s + I/(e*V_a)
## The gain G_p = xi*a((N-N_g)


N_tol = 1
S_tol = 1
for idx in range(n_steps-1):
    dt = dt_next

    #Short Step
    S_short = S_p[idx] + ((G_p[idx] - G_th) * S_p[idx] + xi * a * N[idx]) * dt
    N_short = N[idx] + (-G_p[idx] * S_p[idx] / V_a - N[idx] / tau_s + I / (e * V_a)) * dt

    #Long Step
    S_long = S_p[idx] + ((G_p[idx] - G_th) * S_p[idx] + xi * a * N[idx]) * 2*dt
    N_long = N[idx] + (-G_p[idx] * S_p[idx] / V_a - N[idx] / tau_s + I / (e * V_a)) * 2*dt

    if (np.abs(N_long-N_short) >N_tol) or (np.abs(S_long-S_short) > S_tol):
        dt_next = 0.5 * dt
        N[idx + 1]= N_short
        S_p[idx + 1] = S_short
        step_size[idx + 1] = 0.5*dt
    else:
        dt_next = 2*dt
        N[idx + 1]= N_long
        S_p[idx + 1] = S_long
        step_size[idx + 1] = 2*dt
    G_p[idx + 1] = xi*a*(N[idx] - N_g)
    time[idx + 1] = time[idx] + dt

    #Reset Tolerance based on Value
    N_tol = 0.00001*N[idx]
    S_tol = 0.00001*S_p[idx]

plt.subplot(3, 1, 1)
plt.plot(time,S_p*(hbar*omega/(2*np.pi))*Rb)
plt.subplot(3, 1, 2)
plt.plot(time, N)
plt.plot(N_th*np.ones(np.size(N)))
plt.subplot(3, 1, 3)
plt.plot(step_size)
plt.show()
