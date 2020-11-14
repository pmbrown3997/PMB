function MHD_model()
%% Andrew Knisely
close all; clear; clc;
%%
%% Acceleration of Gravity:
g = 0.001; 
%% Number of Samples on the finite domain:
N = 75;
%% Region of the density:
L = pi/2;  dx = 2*L/N; N = N + 1; dy = dx;
x = -L:dx:L;  y = x; 
X = repmat(x,N,1); Y = repmat(y',1,N);
%%
ro = 4e06;%[per cm^3]
%% Density Region Values:
ro1 = 5*ro;%% Top Region
ro2 = 1*ro;%% Bottom Region 
%%
%% Initial Density Distribution:
%% --------------------------------------------------------------------- %%
RO = zeros(N,N);
for ii = 1:N
    for jj = 1:N
        if(ii < (N+1)/2)
            RO(ii,jj) = ro2;
        else
            RO(ii,jj) = ro1;
        end
    end
end
%% --------------------------------------------------------------------- %%
%%
TIME = 0; 
%% 
%% Initial Perturbation at the Interface:
PSI = .1*sin(X*2+pi).*exp(-4*abs(Y));
%%
%% Calculate the initial Velocity Gradients induced by the perturbation:
Vy = .5*([PSI(:,2:N),PSI(:,1)]-[PSI(:,N),PSI(:,1:N-1)])/dx;   
Vx = -([PSI(2:N,:);0*(1:N)]-PSI)/dy;
%%
%% Gradient for Pressure balanced with the gravitational force:
P = zeros(N,1);
for jj = 1:N
    P(:,jj) = trap_int(RO(:,jj),N,dx)*g;
    P(:,jj)=P(:,jj)-P(N,jj);
end;
%%
dt = 0.16;%% Time Steps
%% Calculate the Initial density distribution from the intial perturbation
RO = advection(1,X,Y,RO,Vx,Vy,dx,dt,N); 
%%
for cnt = 0:20000
   %%
   figure(1); surf(X,Y,RO); shading interp; view(2); axis([-L L -L L]);
   TIME=TIME+dt; str = sprintf('Rayleigh Taylor Instability @ t=%5.1f',TIME);
   colorbar; colormap('jet'); title(str); xlabel('X'); ylabel('Y');
    %% ----------------------------------------------------------------- %%
    %% Calculation of the velocities using momentum eqns:
    %% ----------------------------------------------------------------- %%
    Vx = Vx - ((1./RO).*(.5*([P(:,2:N),P(:,1)]-...
                                              [P(:,N),P(:,1:N-1)])/dx))*dt;
    Vy = Vy-g*dt-((1./RO).*([P(2:N,:);0*(1:N)]-P)/dy)*dt;
    %% ----------------------------------------------------------------- %%
    %% Apply Hodge Decomposition:
    [Vx,Vy] = HHD(Vx,Vy,N,dx,L);
    %%
    %% Update the Density and Velocity Vectors, boundary conditions:
    RO = advection(1,X,Y,RO,Vx,Vy,dx,dt,N); 
    Vx = advection(2,X,Y,RO,Vx,Vy,dx,dt,N);  
    Vy = advection(3,X,Y,RO,Vx,Vy,dx,dt,N);
    %% Boundary Conditions:
    Vy(N,:)=0;         Vy(1,:)=0; 
    %% Update the Gradient for Pressure balanced with the gravitational force:
    for qq = 1:N
        P(:,qq) = trap_int(RO(:,qq),N,dx)*g;
        P(:,qq) = P(:,qq)-P(N,qq);
    end;
end;
%% Trapezoidal Integration Function:
function V = trap_int(a,N,dt)
f = a;   V = zeros(1,N); x = linspace(0,N*dt,N);
   for kk = 1:N
     V(kk) = trapz(x(kk),f(kk),1);
   end
end
%%
end