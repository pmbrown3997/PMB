%% 1D MPD Flow Scalar Conductivity

clear all
clc

% Cylindrical MPD Thruster with Purely radial Current:
% J = j*2pir*z0
%1D steady state plasma acceleration
%% Initial and Constant Values

%Boltzmann Constant
k = 1.38*10^(-23); 

% Number Density of Xenon
n = 10*10^20;%cm^-3
%n = 1000000000000000;%m^-3

%Mass of Xenon
m = 2.18*10^-25; %kg

%vaccuum permiability
mu0 = 1.256*10^(-6);

%Length of Conical Cathode
z_cathode_tip = 0.2; %m lengnth of Cathode tip

%MPD Length
z_max = z_cathode_tip; %m lengnth of MPD thruster

x_max = 0.051;
y_max = 0.051;

% Maximum Radius of Cylindrical MPD
r_max = sqrt(x_max^2 + y_max^2);

%A = pi*r_max^2; %Cross Sectional Area of MPD Thruster

%rho0 = rho0_T*A;%at z=0
%sigma = 55000;

%Initial Density 
rho0 = m*n;

% Gamma = 5/3
gamma = 1.667;

%Inital Preassure
p0 = 10^4; %Pa

%Current
J = 1500; %Amps

% Initial adiabatic sound speed
a0 = sqrt(gamma*(p0/rho0));

% Initial Velocity
u0 = 2000; %m/s

%Initial Mach Number
M0 = u0/a0;

%Ma0 = M0*a0;

% rho*u = F = constant in 1D steady state
F = u0 * rho0;

% Res = 1/sigma;
% 
% V = J*Res;
% 
% E = V/r_max;

% x and y positions are constant 
x = 0.025; %m
y = 0.025; %m

r = sqrt(x^2 + y^2);

% Initial Temperature
T0 = (m/k)*(p0/rho0);

%Permiativitty of Free Space
eps0 = 8.85*10^(-12);

% Charge of Electron
e = 1.6*10^(-19);


%% Induce Magnetic Field from Applied Current in MPD

% Magnetic Field Induced by pureley radial current density from conical
% cathode
B_theta = @(J,r,z_cathode_tip,z) ((mu0 * J)/(2*pi*r))*(1-(z/z_cathode_tip));

%Spitzer-Harm Scalar Conductivity
sigma = @(T) (0.0153*T^(3/2))/(log((12*sqrt(2*pi)*(k*eps0*T)^(3/2))/(e^(3) * sqrt(n))));

%Applied Constant Electric Field
E = @(T) 100; % Alternative Electric Field that treats plasma as basic resistor % (1/r_max)*J*(1/sigma(T));

%% Solution matices

% Length of MPD
L = z_max; %m

deltaz = 0.00001; %Step size in z

N = L/deltaz; %Number of Steps


z = 0:deltaz:L-deltaz; %z array

% Initial Conductivity 
sigma0 = sigma(T0);

u = zeros(1,N);
p = zeros(1,N);
M = zeros(1,N);
a = zeros(1,N);
rho = zeros(1,N);
T = zeros(1,N);
Sigma = zeros(1,N);
%Ma = zeros(1,N);

u(1,1) = u0;
p(1,1) = p0;
M(1,1) = M0;
a(1,1) = a0;
rho(1,1) = rho0;
T(1,1) = T0;
Sigma(1,1) = sigma0;
%Ma(1,1) = Ma0;

f1 = (gamma-1)/gamma;

%% Differential Equations

%Differential Equations obtained for 1D steady state MHD plasma
%acceleration using ideal gas equation of state: 
%c_p * T = (gamma/(gamma-1)) * p/rho

du_dz = @(z,u,p,T,a) ((-mu0*sigma(T)*(f1*E(T)/(mu0*u) - B_theta(J,r,z_cathode_tip,z)/mu0))/(F*(1-(f1 + p/(u*F))))*(E(T)-u*B_theta(J,r,z_cathode_tip,z)));


dp_dz = @(z,u,p,T,a) -(f1*F + p/u)*((-mu0*sigma(T)*(f1*E(T)/(mu0*u) - B_theta(J,r,z_cathode_tip,z)/mu0))/(F*(1-(f1 + p/(u*F))))*(E(T)-u*B_theta(J,r,z_cathode_tip,z)))...
                    - f1*E(T)/(mu0*u)*(-mu0*sigma(T)*(E(T)-u*B_theta(J,r,z_cathode_tip,z)));
                
                
dT_dz = @(z,u,p,T,a) (m/(k*F))*(u*(-(f1*F + p/u)*((-mu0*sigma(T)*(f1*E(T)/(mu0*u) - B_theta(J,r,z_cathode_tip,z)/mu0))/(F*(1-(f1 + p/(u*F))))*(E(T)-u*B_theta(J,r,z_cathode_tip,z)))...
                    - f1*E(T)/(mu0*u)*(-mu0*sigma(T)*(E(T)-u*B_theta(J,r,z_cathode_tip,z))))...
                            + p*(((-mu0*sigma(T)*(f1*E(T)/(mu0*u) - B_theta(J,r,z_cathode_tip,z)/mu0))/(F*(1-(f1 + p/(u*F))))*(E(T)-u*B_theta(J,r,z_cathode_tip,z)))));


da_dz = @(z,u,p,T,a) (gamma/(2*a*F))*(u*(-(f1*F + p/u)*((-mu0*sigma(T)*(f1*E(T)/(mu0*u) - B_theta(J,r,z_cathode_tip,z)/mu0))/(F*(1-(f1 + p/(u*F))))*(E(T)-u*B_theta(J,r,z_cathode_tip,z)))...
                    - f1*E(T)/(mu0*u)*(-mu0*sigma(T)*(E(T)-u*B_theta(J,r,z_cathode_tip,z))))...
                            + p*(((-mu0*sigma(T)*(f1*E(T)/(mu0*u) - B_theta(J,r,z_cathode_tip,z)/mu0))/(F*(1-(f1 + p/(u*F))))*(E(T)-u*B_theta(J,r,z_cathode_tip,z)))));


%%

%RK4 Method for Solving ODE's

% Sonic Point Index
X = 0;

for i = 1:N-1
    
    k1 = du_dz(z(1,i),u(1,i),p(1,i),T(1,i),a(1,i)); 
    w1 = dp_dz(z(1,i),u(1,i),p(1,i),T(1,i),a(1,i)); 
    h1 = dT_dz(z(1,i),u(1,i),p(1,i),T(1,i),a(1,i)); 
    e1 = da_dz(z(1,i),u(1,i),p(1,i),T(1,i),a(1,i)); 

    k2 = du_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*k1*deltaz,p(1,i)+0.5*k1*deltaz,T(1,i)+0.5*k1*deltaz,a(1,i)+0.5*k1*deltaz); 
    w2 = dp_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*w1*deltaz,p(1,i)+0.5*w1*deltaz,T(1,i)+0.5*w1*deltaz,a(1,i)+0.5*w1*deltaz); 
    h2 = dT_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*h1*deltaz,p(1,i)+0.5*h1*deltaz,T(1,i)+0.5*h1*deltaz,a(1,i)+0.5*h1*deltaz); 
    e2 = da_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*e1*deltaz,p(1,i)+0.5*e1*deltaz,T(1,i)+0.5*e1*deltaz,a(1,i)+0.5*e1*deltaz); 

    
    k3 = du_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*k2*deltaz,p(1,i)+0.5*k2*deltaz,T(1,i)+0.5*k2*deltaz,a(1,i)+0.5*k2*deltaz); 
    w3 = dp_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*w2*deltaz,p(1,i)+0.5*w2*deltaz,T(1,i)+0.5*w2*deltaz,a(1,i)+0.5*w2*deltaz);
    h3 = dT_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*h2*deltaz,p(1,i)+0.5*h2*deltaz,T(1,i)+0.5*h2*deltaz,a(1,i)+0.5*h2*deltaz); 
    e3 = da_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*e2*deltaz,p(1,i)+0.5*e2*deltaz,T(1,i)+0.5*e2*deltaz,a(1,i)+0.5*e2*deltaz); 

    
    k4 = du_dz(z(1,i)+0.5*deltaz,u(1,i)+k3*deltaz,p(1,i)+k3*deltaz,T(1,i)+k3*deltaz,a(1,i)+k3*deltaz); 
    w4 = dp_dz(z(1,i)+0.5*deltaz,u(1,i)+w3*deltaz,p(1,i)+w3*deltaz,T(1,i)+w3*deltaz,a(1,i)+w3*deltaz); 
    h4 = dT_dz(z(1,i)+0.5*deltaz,u(1,i)+h3*deltaz,p(1,i)+h3*deltaz,T(1,i)+h3*deltaz,a(1,i)+h3*deltaz); 
    e4 = da_dz(z(1,i)+0.5*deltaz,u(1,i)+e3*deltaz,p(1,i)+e3*deltaz,T(1,i)+e3*deltaz,a(1,i)+e3*deltaz); 


    
    u(1,i+1) = u(1,i) + deltaz*(k1 + 2*k2 + 2*k3 + k4)/6;
    p(1,i+1) = p(1,i) + deltaz*(w1 + 2*w2 + 2*w3 + w4)/6;
    T(1,i+1) = T(1,i) + deltaz*(h1 + 2*h2 + 2*h3 + h4)/6;
    a(1,i+1) = a(1,i) + deltaz*(e1 + 2*e2 + 2*e3 + e4)/6;


    
    rho(1,i+1) = F/u(1,i+1);
    %a(1,i+1) = sqrt(gamma*p(1,i+1)/rho(1,i+1));
    M(1,i+1) = u(1,i+1)/a(1,i+1);
    %T(1,i+1) = (m/k)*(p(1,i+1)/rho(1,i+1));
    Sigma(1,i+1) = sigma(T(1,i+1));
    
    % Sonic Point Check
    if M(1,i+1) >= 1
        X = i+1;
        fprintf('Sonic Point Reached')
        break
    else
        X = 0;
        
    end
    
end

%% Crossing Sonic Point

if X>0
    % subscript 1 means sub-sonic side and 2 means supersonic side
    
    M1 = M(1,X-1);
    u1 = u(1,X-1);
    a1 = a(1,X-1);
    T1 = T(1,X-1);
    z1 = z(1,X-1);
    p1 = p(1,X-1);
    
    
    dudz1 = du_dz(z1,u1,p1,T1,a1);
    dTdz1 = dT_dz(z1,u1,p1,T1,a1);
    dpdz1 = dp_dz(z1,u1,p1,T1,a1);
    dadz1 = da_dz(z1,u1,p1,T1,a1);
    
    M2 = 2-M1;
    
    f = (10/3)*(k/m);
    %Initial valus for supersonic side
    
    T2 = T1;
    
    
    threshold = 0.001;
    
    
    iteration_limit = threshold*200000;
    
    for i = 1:iteration_limit
        
        T2_prev = T2;
        
        
        u2 = M2*(f*(T2_prev))^(1/2);
        
        delz = (u2 - u1)/dudz1;
        
        T2 = T1 + dTdz1*delz;
        
        p2 = p1 + dpdz1*delz;
        
        a2 = a1 + dadz1*delz;
        %     rho2 = F/u2;
        %     a2 = u2/M2;
        %
        %     p2 = (a2^2)*(rho2/gamma);
        %
        z2 = z1 + delz;
        
        q = abs(T2 - T2_prev);
        
        iteration_limit = iteration_limit + 1;
        
        if q < 0.001
            z(1,X) = z2;
            T(1,X) = T2;
            u(1,X) = u2;
            p(1,X) = p2;
            a(1,X) = a2;
            M(1,X) = M2;
            rho(1,X) = F/u2;
            Sigma(1,X) = sigma(T2);
            break
            %     elseif abs(T2 - T2_prev) > 0.001
            %         iteration_limit = iteration_limit + 1;
        end
        
    end
    
    for i = X:N-1
        
        k1 = du_dz(z(1,i),u(1,i),p(1,i),T(1,i),a(1,i));
        w1 = dp_dz(z(1,i),u(1,i),p(1,i),T(1,i),a(1,i));
        h1 = dT_dz(z(1,i),u(1,i),p(1,i),T(1,i),a(1,i));
        e1 = da_dz(z(1,i),u(1,i),p(1,i),T(1,i),a(1,i));
        
        k2 = du_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*k1*deltaz,p(1,i)+0.5*k1*deltaz,T(1,i)+0.5*k1*deltaz,a(1,i)+0.5*k1*deltaz);
        w2 = dp_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*w1*deltaz,p(1,i)+0.5*w1*deltaz,T(1,i)+0.5*w1*deltaz,a(1,i)+0.5*w1*deltaz);
        h2 = dT_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*h1*deltaz,p(1,i)+0.5*h1*deltaz,T(1,i)+0.5*h1*deltaz,a(1,i)+0.5*h1*deltaz);
        e2 = da_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*e1*deltaz,p(1,i)+0.5*e1*deltaz,T(1,i)+0.5*e1*deltaz,a(1,i)+0.5*e1*deltaz);
        
        
        k3 = du_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*k2*deltaz,p(1,i)+0.5*k2*deltaz,T(1,i)+0.5*k2*deltaz,a(1,i)+0.5*k2*deltaz);
        w3 = dp_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*w2*deltaz,p(1,i)+0.5*w2*deltaz,T(1,i)+0.5*w2*deltaz,a(1,i)+0.5*w2*deltaz);
        h3 = dT_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*h2*deltaz,p(1,i)+0.5*h2*deltaz,T(1,i)+0.5*h2*deltaz,a(1,i)+0.5*h2*deltaz);
        e3 = da_dz(z(1,i)+0.5*deltaz,u(1,i)+0.5*e2*deltaz,p(1,i)+0.5*e2*deltaz,T(1,i)+0.5*e2*deltaz,a(1,i)+0.5*e2*deltaz);
        
        
        k4 = du_dz(z(1,i)+0.5*deltaz,u(1,i)+k3*deltaz,p(1,i)+k3*deltaz,T(1,i)+k3*deltaz,a(1,i)+k3*deltaz);
        w4 = dp_dz(z(1,i)+0.5*deltaz,u(1,i)+w3*deltaz,p(1,i)+w3*deltaz,T(1,i)+w3*deltaz,a(1,i)+w3*deltaz);
        h4 = dT_dz(z(1,i)+0.5*deltaz,u(1,i)+h3*deltaz,p(1,i)+h3*deltaz,T(1,i)+h3*deltaz,a(1,i)+h3*deltaz);
        e4 = da_dz(z(1,i)+0.5*deltaz,u(1,i)+e3*deltaz,p(1,i)+e3*deltaz,T(1,i)+e3*deltaz,a(1,i)+e3*deltaz);
        
        
        
        u(1,i+1) = u(1,i) + deltaz*(k1 + 2*k2 + 2*k3 + k4)/6;
        p(1,i+1) = p(1,i) + deltaz*(w1 + 2*w2 + 2*w3 + w4)/6;
        T(1,i+1) = T(1,i) + deltaz*(h1 + 2*h2 + 2*h3 + h4)/6;
        a(1,i+1) = a(1,i) + deltaz*(e1 + 2*e2 + 2*e3 + e4)/6;
        
        
        
        rho(1,i+1) = F/u(1,i+1);
        %a(1,i+1) = sqrt(gamma*p(1,i+1)/rho(1,i+1));
        M(1,i+1) = u(1,i+1)/a(1,i+1);
        %T(1,i+1) = (m/k)*(p(1,i+1)/rho(1,i+1));
        Sigma(1,i+1) = sigma(T(1,i+1));
        
        if M(1,i) - M(1,i-1) >1
            break
        end
        
    end
    
    
end


%% PLot Solutions



figure(1)
subplot(2,2,1)
plot(z,u)
grid on
xlabel('ditance along MPD (m)')
ylabel('velocity (m/s)')
title('Velocity Along MPD')

subplot(2,2,2)
plot(z,M)
grid on
xlabel('ditance along MPD (m)')
ylabel('Mach number')
title('Mach Number Along MPD')

subplot(2,2,3)
plot(z,rho)
grid on
xlabel('ditance along MPD (m)')
ylabel('linear density (kg/l)')
title('Linear Density Along MPD')

subplot(2,2,4)
plot(z,p)
grid on
xlabel('ditance along MPD (m)')
ylabel('linear preassure (Pa)')
title('Linear Preassure Along MPD')


figure(2)
plot(z,Sigma)
grid on
xlabel('ditance along MPD (m)')
ylabel('Conductivity (SI Units)')
title('Conductivity Along MPD')
