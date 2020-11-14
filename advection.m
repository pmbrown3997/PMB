function T=advection(type,X,Y,RO,Vx,Vy,dx,dt,N)
%% Periodic BC in X and fixed BC in Y
%%
%% Calculate the advection of each quantity:
%type refers to 1. Density, 2. Vx, 3. Vy.
%%
if(type==1) 
    T=RO; 
elseif(type==2) 
    T=Vx; 
else
    T=Vy; 
end;
%%
T0=T;   xmin=X(1,1);      ymin=Y(1,1); 
%%
%% Calculate the positions from initial velocity:
dX=Vx*dt;                                 dY=Vy*dt;
Xold=X-xmin-dX;                           Yold = Y-ymin-dY;
%%
%% Index:
Iold=floor(Xold/dx)+1;                    Jold=floor(Yold/dx)+1;
fx = Xold/dx-Iold+1;                      fy = Yold/dx-Jold+1;
F1 = (1-fx).*(1-fy);                      F2 = fx.*(1-fy);
F3 =  fx.*fy;                             F4 = (1-fx).*fy;
Iold=Iold.*(Iold<=N)+N.*(Iold>N);         Jold=Jold.*(Jold<=N)+N.*(Jold>N);
%%
for ix=1:N
    for iy = 1:N
        %%
        ix1 = Iold(iy,ix); iy1 = Jold(iy,ix);
        %%
        if(ix1>0&&ix1<N&&iy1>0&&iy1<N)
            ro1 = T0(iy1,ix1);       ro2 = T0(iy1,ix1+1);
            ro3 = T0(iy1+1,ix1+1);   ro4 = T0(iy1+1,ix1);
        end;
        %%
        if(ix1<=0)
            if(iy1<N&&iy1>0)
               ro1 = T0(iy1,N);     ro2 = T0(iy1,1);
               ro3 = T0(iy1+1,1);   ro4 = T0(iy1+1,N);
            end;
            %%
            if(iy1<=0)
                ro3 = T0(1,1);      ro4 = T0(1,N);
                ro1 = ro4;          ro2 = ro3;
            end
            %%
            if(iy1>=N)
                ro1 = T0(N,N);      ro2 = T0(N,1);
                ro3 = ro2;          ro4 = ro1;
             end;
        end;
        %%
        if(ix1>=N)
            if(iy1<N&&iy1>0)
               ro1 = T0(iy1,N);    ro2 = T0(iy1,1);
               ro3 = T0(iy1+1,1);  ro4 = T0(iy1+1,N);
            end;
        %%
            if(iy1<=0)
                ro3 = T0(1,1);     ro4 = T0(1,N);
                ro1 = ro4;         ro2 = ro3;
            end
        %%
            if(iy1>=N)
                ro1 = T0(N,N);     ro2 = T0(N,1);
                ro3 = ro2;         ro4 = ro1;
             end;
        %%
        end;
        if(iy1<=0&&ix1>0&&ix1<N)
            ro3 = T0(1,ix1);      ro4 = T0(1,ix1+1);
            ro1 = ro4;            ro2 = ro3;
        end;
        %%
        if(iy1>=N&&ix1>0&&ix1<N)
            ro1 = T0(N,ix1);      ro2 = T0(N,ix1+1);
            ro4 = ro1;            ro3 = ro2;
        end;
        T(iy,ix) = ro1*F1(iy,ix)+ro2*F2(iy,ix)+...
                                               ro3*F3(iy,ix)+ro4*F4(iy,ix);
        %%
    end;
end;