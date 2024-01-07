clear all;
xi = 0.2;
V_a = 1e-16;
e = 1.6e-19;
c=2.998e8;
n_eq = 3.1;
a = 2.75e-12;
L = 0.001;
Rb = 0.99;
Rf = 0.99;
alpha_loss = 1/8.68;
nu = c/n_eq;
N_g = 2.1e24;
T = 298.2;
n_current = 60;

hbar = 1.055e-34;
omega = 2*pi*3.527e14;
npts = 2^24;

G_th = nu*(alpha_loss +1/(2*L)*log(1/(Rb*Rf)));
N_th = G_th/(xi*a) + N_g;
test = xi*a*(N_th-N_g);
dt = 1e-15;
S_p = eps*ones(1,npts);
N = N_th*ones(1,npts);
%N = zeros(1,npts);
G_p = zeros(1,npts);
N_c = 3.65e18;
Pwr = zeros(npts,n_current);
tau_s = 1e3*2.85e7*T/(N_c*1e6);

I_th = e*N_th*V_a/tau_s;
I = linspace(2*I_th,2*I_th,n_current);
 
S_tol = 100;
N_tol = N_g/1000;

y0 = [0 N_th];
tspan = [0,1e-05];

[t,y] = ode45(@odefun, tspan, y0);
Pwr = nu/(2*L)*(1/(Rb*Rf))*y(:,1)*hbar*omega;
figure()
plot(t,Pwr)
figure()
plot(t,y(:,2))

function dydt = odefun(t,y)
xi = 0.2;
V_a = 1e-16;
e = 1.6e-19;
c=2.998e8;
n_eq = 3.1;
a = 2.75e-12;
L = 0.001;
Rb = 0.99;
Rf = 0.99;
alpha_loss = 1/8.68;
nu = c/n_eq;
N_g = 2.1e24;
T = 298.2;
G_th = nu*(alpha_loss +1/(2*L)*log(1/(Rb*Rf)));
N_th = G_th/(xi*a) + N_g;
%N = zeros(1,npts);
N_c = 3.65e18;
tau_s = 1e3*2.85e7*T/(N_c*1e6);

I = 0.9*e*N_th*V_a/tau_s;

dydt = zeros(2,1);
dydt(1) = (xi*a*(y(2) - N_g) -G_th)*y(1) + xi*a*y(2);
dydt(2) = -xi*a*(y(2) - N_g)*y(1)/V_a -y(2)/tau_s + I/(e*V_a);
end


