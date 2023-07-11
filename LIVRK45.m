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
npts = 2^25;

G_th = nu*(alpha_loss +1/(2*L)*log(1/(Rb*Rf)));
N_th = G_th/(xi*a) + N_g;
test = xi*a*(N_th-N_g);
dt = 1e-15;
t = zeros(1,npts);
S_p = eps*ones(1,npts);
N = N_th*ones(1,npts);
%N = zeros(1,npts);
G_p = zeros(1,npts);
N_c = 3.65e18;
Pwr = zeros(npts,n_current);
tau_s = 1e3*2.85e7*T/(N_c*1e6);

I_th = e*N_th*V_a/tau_s;
I = linspace(2*I_th,2*I_th,n_current);
 
S_tol = 10;
N_tol = N_g/1000;

for current = 1:n_current
    
    for idx = 1:(npts-1)
        %tau_s = 2.85e7*T/(N_c*1e6);
        G_p(idx+1) = xi*a*(N(idx+1) - N_g);

        %First Calc
        N_short = N(idx) + (-G_p(idx)*S_p(idx)/V_a -N(idx)/tau_s + I(current)/(e*V_a))*dt;
        S_short = S_p(idx) +((G_p(idx) -G_th)*S_p(idx) + xi*a*N(idx))*dt;

        %Second Calc 
        N_long = N(idx) + (-G_p(idx)*S_p(idx)/V_a -N(idx)/tau_s + I(current)/(e*V_a))*(2*dt);
        S_long = S_p(idx) +((G_p(idx) -G_th)*S_p(idx) + xi*a*N(idx))*(2*dt);

        if (abs(N_long-N_short)>N_tol) || (abs(S_long-S_short) > S_tol)
            dt = 0.5*dt;
            N(idx+1) = N_short;
            S_p(idx+1) = S_short;
        else
            dt = 2*dt;
            N(idx+1) = N_short;
            S_p(idx+1) = S_short;
        end
        
        S_tol = 0.001*S_p(idx+1);
        N_tol = 0.001*N(idx+1);
        
        t(idx+1) = t(idx)+dt;
        

        %S_p(idx+1) = S_p(idx) + ((G_p(idx) -G_th)*S_p(idx))*dt;
        if(S_p(idx+1)<0)
            S_p(idx+1) = 0;
        end
        
    end
    figure()
    plot(t,S_p)
    figure()
    plot(t, N)
    hold on;
    plot(t,ones(size(N))*N_th);
    Pwr(:,current) = nu/(2*L)*(1/(Rb*Rf))*S_p*hbar*omega; 
end

plot(I,Pwr(end,:));
