clear all
clf
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% FLOW OVER A STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Mar¨ªa Fern¨¢ndez Jim¨¦nez %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EQUATIONS
% Continuity
% p(n+1)=p(n)-(dt/(2*dx))*(u(i+1)-u(i-1))-(dt/(2*dy))*(v(j+1)-v(j-1));
% X-Momentum
% u(n+1)=u(n)-((u*dt)/(2dx))*(u(i+1)-u(i-1))-((v*dt)/(2dy))*(u(j+1)-u(j-1))-(dt*(2*dx))*(p(i+1)-p(i-1))+mu*dt/Re*((u(i+1)+u(i-1)-2*u)/(dx^2)+(u(j+1)+u(j-1)-2*u)/(dy^2));
% Y-Momentum
% v(n+1)=v(n)-((u*dt)/(2dx))*(v(i+1)-v(i-1))-((v*dt)/(2dy))*(v(j+1)-v(j-1))-(dt*(2*dy))*(p(j+1)-p(j-1))+mu*dt/Re*((v(i+1)+v(i-1)-2*v)/(dx^2)+(v(j+1)+v(j-1)-2*v)/(dy^2));
% CFL = max((dt*u/dx , dt*v/dy)) < 1
% COUNTERS
dy=0.01; y=[0:dy:1]; Ny=length(y);
dx=0.01; x=[0:dx:10]; Nx=length(x);
dt=0.0001; Nt=10000000;
% CONSTANTS
Re=100; mu=1;
H=round(Ny/2);
L=round(Nx/3.5);
% MATRICES INITIATION
u_new=zeros(Nx,Ny);
v_new=zeros(Nx,Ny);
p_new=zeros(Nx,Ny);
phi=zeros(Nx,Ny);
contador=1;


for H=round(Ny/2):-10:1
%% INITIAL CONDITIONS
a = y(Ny)-y(H);
bprime = y(round((Ny-H)/2)+H)-y(H);
b = 1.5/(bprime*(bprime-a));
i=1;
for j=1:Ny
u_old(i,j) = b*(y(j)-y(H))*((y(j)-y(H))-a);
v_old(i,j) = 0;
p_old(i,j) = 0;
end
%% PROBLEM RESOLUTION
T=1;
max_error=1e-3;
error=10*max_error;
while error>max_error
for i=2:Nx-1
for j=2:Ny-1
if ((i<L)&&(j<H))
u_new(i,j)=0; v_new(i,j)=0; p_new(i,j)=0;
else
% MOMENTUM & CONTINUITY EQUATIONS
dudx = (u_old(i+1,j)-u_old(i-1,j))/(2*dx);
dvdx = (v_old(i+1,j)-v_old(i-1,j))/(2*dx);
dpdx = (p_old(i+1,j)-p_old(i-1,j))/(2*dx);
dudy = (u_old(i,j+1)-u_old(i,j-1))/(2*dy);
dvdy = (v_old(i,j+1)-v_old(i,j-1))/(2*dy);
dpdy = (p_old(i,j+1)-p_old(i,j-1))/(2*dy);
d2udx2=(u_old(i+1,j)+u_old(i-1,j)- 2*u_old(i,j))/(dx^2);
d2vdx2=(v_old(i+1,j)+v_old(i-1,j)- 2*v_old(i,j))/(dx^2);
d2udy2=(u_old(i,j+1)+u_old(i,j-1)- 2*u_old(i,j))/(dy^2);
d2vdy2=(v_old(i,j+1)+v_old(i,j-1)- 2*v_old(i,j))/(dy^2);
u_new(i,j)=u_old(i,j)-u_old(i,j)*dt*dudxv_old(i,j)*dt*dudy-dt*dpdx+(mu/Re)*dt*(d2udx2+d2udy2);
v_new(i,j)=v_old(i,j)-u_old(i,j)*dt*dvdxv_old(i,j)*dt*dvdy-dt*dpdy+(mu/Re)*dt*(d2vdx2+d2vdy2);
p_new(i,j) = p_old(i,j) - dt*dudx - dt*dvdy;
end
end
end
%% BOUNDARY CONDITIONS
% Upper Wall
for i=2:Nx-1
j=Ny;
u_new(i,j) = 0; 
v_new(i,j) = 0; 
d2vdx2 = (v_new(i+1,j)+v_new(i-1,j)-2*v_new(i,j)) / (dx^2);
d2vdy2_up = (v_new(i,j)+v_new(i,j-2)-2*v_new(i,j-1)) / (dy^2);
p_new(i,j) = p_new(i,j-1) + (mu/Re)*dy*(d2vdx2 + d2vdy2_up);
end
j=Ny;
i=1;
u_new(i,j) = 0; 
v_new(i,j) = 0; 
d2vdx2_left = (v_new(i,j)+v_new(i+2,j)-2*v_new(i+1,j)) / (dx^2);
d2vdy2_up = (v_new(i,j)+v_new(i,j-2)-2*v_new(i,j-1)) / (dy^2);
p_new(i,j) = p_new(i,j-1) + (mu/Re)*dy*(d2vdx2_left + d2vdy2_up);
i=Nx;
u_new(i,j) = 0; 
v_new(i,j) = 0; 
d2vdx2_right = (v_new(i,j)+v_new(i-2,j)-2*v_new(i-1,j)) / (dx^2);
d2vdy2_up = (v_new(i,j)+v_new(i,j-2)-2*v_new(i,j-1)) / (dy^2);
p_new(i,j) = p_new(i,j-1) + (mu/Re)*dy*(d2vdx2_right + d2vdy2_up);
% Lower Wall After Step
for i=L+1:Nx-1
j=1;
u_new(i,j) = 0; 
v_new(i,j) = 0; 
d2vdx2 = (v_new(i+1,j)+v_new(i-1,j)-2*v_new(i,j)) / (dx^2);
d2vdy2_low = (v_new(i,j)+v_new(i,j+2)-2*v_new(i,j+1)) / (dy^2);
p_new(i,j) = p_new(i,j+1) - (mu/Re)*dy*(d2vdx2 + d2vdy2_low);
end
j=1;
i=L;
u_new(i,j) = 0; 
v_new(i,j) = 0; 
d2vdx2_left = (v_new(i,j)+v_new(i+2,j)-2*v_new(i+1,j)) / (dx^2);
d2vdy2_low = (v_new(i,j)+v_new(i,j+2)-2*v_new(i,j+1)) / (dy^2);
p_new(i,j) = p_new(i,j+1) - (mu/Re)*dy*(d2vdx2_left + d2vdy2_low); 
i=Nx;
u_new(i,j) = 0; 
v_new(i,j) = 0; 
d2vdx2_right = (v_new(i,j)+v_new(i-2,j)-2*v_new(i-1,j)) / (dx^2);
d2vdy2_low = (v_new(i,j)+v_new(i,j+2)-2*v_new(i,j+1)) / (dy^2);
p_new(i,j) = p_new(i,j+1) - (mu/Re)*dy*(d2vdx2_right + d2vdy2_low); 
% Lower Wall Before Step
for i=2:L
j=H+1;
u_new(i,j) = 0; 
v_new(i,j) = 0; 
d2vdx2 = (v_new(i+1,j)+v_new(i-1,j)-2*v_new(i,j)) / (dx^2);
d2vdy2_low = (v_new(i,j)+v_new(i,j+2)-2*v_new(i,j+1)) / (dy^2);
p_new(i,j) = p_new(i,j+1) - (mu/Re)*dy*(d2vdx2 + d2vdy2_low);
end
j=H;
i=1;
u_new(i,j) = 0; 
v_new(i,j) = 0; 
d2vdx2_left = (v_new(i,j)+v_new(i+2,j)-2*v_new(i+1,j)) / (dx^2);
d2vdy2_low = (v_new(i,j)+v_new(i,j+2)-2*v_new(i,j+1)) / (dy^2);
p_new(i,j) = p_new(i,j+1) - (mu/Re)*dy*(d2vdx2_left + d2vdy2_low);
% Vertical Wall
for j=1:H
i=L;
u_new(i,j) = 0; 
v_new(i,j) = 0; 
d2udx2_left = (u_new(i,j)+u_new(i+2,j)-2*u_new(i+1,j)) / (dx^2);
d2udy2_low = (u_new(i,j)+u_new(i,j+2)-2*u_new(i,j+1)) / (dy^2);
p_new(i,j) = p_new(i+1,j)-(mu/Re)*dx*(d2udx2_left + d2udy2_low);
end
% Inlet
for j=H+1:Ny-1
i=1;
u_new(i,j) = b*(y(j)-y(H))*((y(j)-y(H))-a);
v_new(i,j) = 0; 
dudx_left = (u_new(i+1,j)-u_new(i,j)) / dx;
d2udx2_left = (u_new(i,j)+u_new(i+2,j)-2*u_new(i+1,j)) / (dx^2);
d2udy2 = (u_new(i,j+1)+u_new(i,j-1)-2*u_new(i,j)) / (dy^2);
p_new(i,j) = p_new(i+1,j)-(mu/Re)*dx*(d2udx2_left+d2udy2)-u_new(i,j)*dudx_left;
end
i=1;
j=H;
u_new(i,j) = b*(y(j)-y(H))*((y(j)-y(H))-a);
v_new(i,j) = 0;
d2udx2_left = (u_new(i,j)+u_new(i+2,j)-2*u_new(i+1,j)) / (dx^2);
d2udy2_low = (u_new(i,j)+u_new(i,j+2)-2*u_new(i,j+1)) / (dy^2);
p_new(i,j) = p_new(i+1,j) - (mu/Re)*dx*(d2udx2_left + d2udy2_low); 
j=Ny;
u_new(i,j) = b*(y(j)-y(H))*((y(j)-y(H))-a);
v_new(i,j) = 0; 
d2udx2_left = (u_new(i,j)+u_new(i+2,j)-2*u_new(i+1,j)) / (dx^2);
d2udy2_up = (u_new(i,j)+u_new(i,j-2)-2*u_new(i,j-1)) / (dy^2);
p_new(i,j) = p_new(i+1,j) - (mu/Re)*dx*(d2udx2_left + d2udy2_up); 
%Outlet
for j=2:Ny-1
i=Nx;
u_new(i,j) = u_new(i-1,j); 
v_new(i,j) = v_new(i-1,j); 
d2udx2_right = (u_new(i,j)+u_new(i-2,j)-2*u_new(i-1,j))/(dx^2);
d2udy2 = (u_new(i,j+1)+u_new(i,j-1)-2*u_new(i,j)) / (dy^2);
p_new(i,j) = p_new(i-1,j) + (mu/Re)*dx*(d2udx2_right + d2udy2);
end
i=Nx;
j=1;
u_new(i,j) = u_new(i-1,j); 
v_new(i,j) = v_new(i-1,j);
d2udx2_right = (u_new(i,j)+u_new(i-2,j)-2*u_new(i-1,j)) / (dx^2);
d2udy2_low = (u_new(i,j)+u_new(i,j+2)-2*u_new(i,j+1)) / (dy^2);
p_new(i,j) = p_new(i-1,j) + (mu/Re)*dx*(d2udx2_right + d2udy2_low);
j=Ny;
u_new(i,j) = u_new(i-1,j); 
v_new(i,j) = v_new(i-1,j);
d2udx2_right = (u_new(i,j)+u_new(i-2,j)-2*u_new(i-1,j)) / (dx^2);
d2udy2_up = (u_new(i,j)+u_new(i,j-2)-2*u_new(i,j-1)) / (dy^2);
p_new(i,j) = p_new(i-1,j) + (mu/Re)*dx*(d2udx2_right + d2udy2_up); 
% CFL Condition Checking
if (T==1)
Rx = dt*max(max(u_new))/dx;
Ry = dt*max(max(v_new))/dy;
CFL = max(Rx,Ry);
if(CFL>=1)
disp(['CFL not valid'])
break
end
end
% Error Calculation for Convergence
errorU = max(max((u_new-u_old)/dt));
errorV = max(max((v_new-v_old)/dt));
errorP = max(max((p_new-p_old)/dt));
mass_inlet = trapz(y,u_new(1,:));
mass_outlet = trapz(y,u_new(end,:));
errormass=abs(mass_inlet-mass_outlet)/mass_inlet;
error=errorU+errorV+errorP+errormass;
% Update Solution
u_old=u_new;
v_old=v_new;
p_old=p_new;
% Print Data in the Screen
if (rem(T,1000)==0)
for i=1:Nx
phi(i,:)=cumtrapz(y,u_new(i,:));
end
% Screen displays
disp(['For Time:'])
disp(T)
disp(['The errors are: (P U V)'])
disp([errorP errorU errorV])
disp(['Mass in the inlet'])
disp(mass_inlet)
disp(['Mass in the outlet'])
disp(mass_outlet)
disp(['Re='])
disp(Re)
disp(['H='])
disp(H)
if (contador==1)
save fileRe100_1.mat u_new v_new p_new x y H L Re
end
if (contador==2)
save fileRe100_2.mat u_new v_new p_new x y H L Re
end
if (contador==3)
save fileRe100_3.mat u_new v_new p_new x y H L Re
end
if (contador==4)
save fileRe100_4.mat u_new v_new p_new x y H L Re
end
if (contador==5)
save fileRe100_5.mat u_new v_new p_new x y H L Re
end
end
T=T+1;
end
contador=contador+1;
end
D=[2:-0.005:-0.1];
for i=1:Nx
phi(i,j)=cumtrapz(y,u_new(i,:));
end
contour(x,y,phi',D)