using Plots

# Define Variables
xi = 0.2                # Field Confinment Factor
V_a = 1e-16             # Volume of Active region
e = 1.6e-19             # Charge of an Electron in Coulombs
c = 2.998e8             # Speed of light in m/s
n_eq = 3.1              # Effective index of lasing mode.
a = 2.75e-12            # Tangential Coefficient
L = 0.001               # Cavity length (meters)
Rb = 0.9                # First Mirror reflectivity  (Intensity)
Rf = 0.9                # Second Mirror reflectivity (Intensity)
alpha_loss = 1/8.68     # Round trip loss coefficient
nu = c/n_eq             # Propagation Speed of optical wave in laser (speed of light / effective index)
N_g = 2.1e24            # Transparent electron density
T = 298.2
N_c = 3.65e18
tau_s = 1e3*2.85e7*T/(N_c*1e6)

hbar = 1.055e-34
omega = 2*pi*3.527e14

t_span = (0, 1e-8)

dt = tau_s/20000

G_th = nu*(alpha_loss + 1/(2*L)*log(1/(Rb*Rf)))
N_th = G_th/(xi*a) + N_g
I_th = e*N_th*V_a/tau_s
I = 1.2*I_th

# Simulation
# The following simulation attempts to solve the following rate equations for Photon number and electron density
# dS_p/dt = (G_p - G_th)S_p + xi*a*N. Note xi*a*N is an approximation of the spontaneous emission term
# dN/dt = -G_p S_p/V_a - N/tau_s + I/(e*V_a)
# The gain G_p = xi*a((N-N_g)

function rk45_step(fun, t, y, h, args...)
    k1 = h * fun(t, y, args...)
    k2 = h * fun(t + 0.2 * h, y + 0.2 * k1, args...)
    k3 = h * fun(t + 0.3 * h, y + 0.075 * k1 + 0.225 * k2, args...)
    k4 = h * fun(t + 0.6 * h, y + 0.3 * k1 + (-0.9) * k2 + 1.2 * k3, args...)
    k5 = h * fun(t + h, y + (-11 / 54) * k1 + 2.5 * k2 + (-70 / 27) * k3 + (35 / 27) * k4, args...)
    k6 = h * fun(t + 0.875 * h, y + (1631 / 55296) * k1 + (175 / 512) * k2 + (575 / 13824) * k3 +
               (44275 / 110592) * k4 + (253 / 4096) * k5, args...)

    y_next = y + (37 / 378) * k1 + (250 / 621) * k3 + (125 / 594) * k4 + (512 / 1771) * k6
    y_star = y + (2825 / 27648) * k1 + (18575 / 48384) * k3 + (13525 / 55296) * k4 + (277 / 14336) * k5 + (1 / 4) * k6

    error = norm(y_next - y_star)
    return y_next, error
end

function rate_equations(t, y, xi, a, N_g, G_th, V_a, tau_s, I, e)
    N, S = y
    dN_dt = -xi*a*(N - N_g) * S / V_a - N / tau_s + I / (e * V_a)
    dS_dt = (xi*a*(N - N_g) - G_th) * S + xi * a * N
    return [dN_dt, dS_dt]
end

function solve_rk45(fun, t_span, y0, h, args...)
    t_values = [t_span[1]]
    y_values = [y0]

    while t_values[end] > t_span[1]
        t = t_values[end]
        y = y_values[end]
        y_next, error = rk45_step(fun, t, y, -h, args...)

        while error > 1e-4
            h /= 2
            y_next, error = rk45_step(fun, t, y, -h, args...)
        end

        push!(t_values, t + h)
        push!(y_values, y_next)
    end

    return reverse(t_values), hcat(reverse(y_values)...)
end

# Initial conditions
N0 = 0.0  # Initial carrier density
S0 = 0.0  # Initial photon number

# Solve the rate equations using RK45 method
t_values, y_values = solve_rk45(rate_equations, t_span, [N0, S0], dt, xi, a, N_g, G_th, V_a, tau_s, I, e)

N = y_values[1, :]
S_p = y_values[2, :]

plot(t_values, S_p .* (hbar*omega/(2*pi)) .* Rb, subplot=2, xlabel="Time", ylabel="Photon Number")
plot!(t_values, N, subplot=1, xlabel="Time", ylabel="Carrier Density")