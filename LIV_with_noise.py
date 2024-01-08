import numpy as np
import matplotlib.pyplot as plt
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

t_span = (0,500e-13)

dt = tau_s/20000000



#G_th_1 = nu*(alpha_loss +1/(2*L)*np.log(1/(Rb*Rf)))
G_th = 5.01e11
N_th = G_th/(xi*a) + N_g
I_th = e*N_th*V_a/tau_s
I = 1.5*I_th

def step(fun, t, y, h,N_avg, *args):
    k1 = h * fun(t, y, h, N_avg, *args)
    y_next = y + k1
    return y_next


def rate_equations(t,y,h,N_avg,xi,a,N_g,G_th,V_a,tau_s,I,e,alpha):
    N = y[0]
    N_g2 = N_g*V_a
    if y[1] <0:
        S = 0
    else:
        S = y[1]
    T = y[2]
    #Calculate V matrix
    V_SS = (2*a*xi/V_a)*N*(S+1)
    V_NN = 2*(1 + a*xi*tau_s/V_a*S)*N/tau_s
    V_SN = -a*xi/V_a*((N + N_g2)*S + N)
    V_TT = V_SS/(4*(S+1)**2)
    V_NT = V_SN/(2*(S + 1))

    V_SS_t = ((a*xi/V_a)*(N+N_g2) + G_th)*S + a*xi*N/V_a
    V_NN_t = a*xi/V_a*(N+N_g2)*S + N/tau_s + I/e
    V_SN_t = -a*xi/V_a*((N+N_g2)*S + N)
    V_TT_t = V_SS_t / (4 * (S + 1) ** 2)
    V_NT_t = V_SN_t / (2 * (S + 1))
    #Calculate k and m
    k = -V_SN/V_SS
    m = -V_NT/V_TT

    k_t = -V_SN_t/V_SS_t
    m_t = -V_NT_t/V_TT_t

    #Estimate Random Increment
    mu, sigma = 0, 1  # mean and standard deviation
    g_s = np.random.normal(mu, sigma, 1)
    g_t = np.random.normal(mu, sigma, 1)
    g_n = np.random.normal(mu, sigma, 1)

    #Calculate Random Increment
    #test = (k**2)*V_SS + k*V_SN + (m**2)*V_TT + m*V_NT + k*V_SN + m*V_NT + V_NN
    #test2 = (k_t**2)*V_SS_t + k_t*V_SN_t + (m_t**2)*V_TT_t + m_t*V_NT_t + k_t*V_SN_t + m_t*V_NT_t + V_NN_t
    F_S = np.sqrt(V_SS/h)*g_s
    F_T = 1/(2*S+1)*np.sqrt(V_SS/h)*g_t


    test_1 = V_NN+2*k*V_SN
    test_2 = -(V_NN_t+2*k_t*V_SN_t)
    F_temp = np.sqrt(np.abs((V_NN+2*k*V_SN))/h)*g_n
    F_N = F_temp - k*F_S - m*F_T

    dN_dt = -xi/V_a*a*(N - N_g2)*S - N/tau_s + I/e + F_N
    dS_dt = (xi/V_a*a*(N - N_g2) - G_th) * S + xi * a * N/V_a + F_S
    #dN_dt = np.array([-xi * a * (N / V_a - N_g2) * S - N / tau_s + I / e])
    #dS_dt = np.array([(xi * a * (N / V_a - N_g2) - G_th) * S + xi * a * N / V_a])

    dT_dt = alpha*a*xi/(2*V_a)*(N-N_avg)+F_T
    return np.array([dN_dt[0], dS_dt[0], dT_dt[0]])



def solve_ode(fun, t_span, y0, h, *args):
    t_values = [t_span[0]]
    y_values = np.array([y0])
    while t_values[-1] < t_span[1]:
        t = t_values[-1]
        y = y_values[-1,:]
        N_avg = np.mean(y_values[:,0])
        y_next = np.array([step(fun, t, y, h, N_avg, *args)])
        t_values.append(t + h)
        y_values = np.vstack([y_values, y_next])
    return np.array(t_values), np.array(y_values).T


# Initial conditions
N0 = N_th*V_a # Initial carrier density
S0 = 4.45990821e9-210+170  # Initial photon number
T0 = 0.0
# Solve the rate equations using RK45 method
t_values, y_values = solve_ode(rate_equations, t_span, np.array([N0, S0, T0]),dt,xi,a,N_g,G_th,V_a,tau_s,I,e,alpha)

N = y_values[0]
S_p = y_values[1]

plt.subplot(2, 1, 1)
plt.plot(t_values,S_p)
plt.subplot(2, 1, 2)
plt.plot(t_values, N)
plt.plot(t_values,N_th*V_a*np.ones(np.size(N)))
plt.show()
