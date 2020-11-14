function [Ux,Uy] = HHD(Vx,Vy,N,dx,L)
%% Helmholtz Decomposition:
    dxVx=([Vx(:,2:N),Vx(:,1)]-[Vx(:,N),Vx(:,1:N-1)])/2/dx;
    dyVy=([Vy(2:N,:);Vy(1,:)]-[Vy(N,:);Vy(1:N-1,:)])/2/dx;
    divV=dxVx+dyVy;
%% 
%% Using the Spectral Method to solve the Poisson Eqn for PHI:
%% Delsq PHI = Div(U):
%%
  N = length(divV);%% No. of Fourier modes...should be a power of 2
  Lr = 2*L;%% Domain size (assumed square)
  k = (2*pi/Lr)*[0:(N/2-1) (-N/2):(-1)];%% Vector of wavenumbers
  [KX,KY]  = meshgrid(k,k); %% Matrix of (x,y) wavenumbers corresponding
                            %% to Fourier mode (m,n)
  delsq = -(KX.^2 + KY.^2); %% Laplacian matrix acting on the wavenumbers
  delsq(1,1) = 1;%% Kluge to avoid division by zero of (0,0) waveno.
%% Construct RHS f(x,y) at the Fourier gridpoints
  f = divV;
%% Spectral inversion of Laplacian
  fhat = fft2(f);
  PHI = real(ifft2(fhat./delsq));
  PHI = PHI - PHI(1,1);%% Specify arbitrary constant by forcing corner u = 0.
%%
  Ux=([PHI(:,2:N),PHI(:,1)]-[PHI(:,N),PHI(:,1:N-1)])/2/dx;
%% Get U=-grad PHI
  Ux=Vx-Ux;
%% Eliminate the compressible flow:
  Uy=([PHI(2:N,:);PHI(1,:)]-[PHI(N,:);PHI(1:N-1,:)])/2/dx;
  Uy(1,:)=0;Uy(N,:)=0; Uy=Vy-Uy;
%%