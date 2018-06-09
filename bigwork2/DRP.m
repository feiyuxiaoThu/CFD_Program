%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=256;
L=1;
a=1;
dx=L/N;
x=0:dx:L;
CFL=0.1;
dt=a*CFL*dx;
t=10*dt;

k0=15;

u=zeros(1,N+1);
uu=zeros(1,N+1);
u(:)=1;
uu(:)=1;
seed=zeros(N+1,N/4);

% for i=1:N+1
%     u(i)=3*sin(2*pi*x(i))+4*sin(100*pi*x(i))+5*sin(500*pi*x(i));
% end
% for i=1:N+1
%     uu(i)=3*sin(2*pi*(x(i)-t))+4*sin(100*pi*(x(i)-t))+5*sin(500*pi*(x(i)-t));
% end
% for i=1:N+1
%     if x(i)<0.25
%         u(i)=0;
%     end
%     if x(i)>0.75
%         u(i)=0;
%     end
% end
for i=1:N+1
    for j=1:N/4
        E=(j/k0)^4*exp(-2*(j/k0)^2);
        seed(i,j)=rand(1,1);
    u(i)=u(i)+0.1*(E^0.5*sin(2*pi*j*(x(i)+seed(i,j))));
    end
end
for i=1:N+1
    for j=1:N/4
        E=(j/k0)^4*exp(-2*(j/k0)^2);
    uu(i)=uu(i)+0.1*(E^0.5*sin(2*pi*j*(x(i)-t+seed(i,j))));
    end
end

figure;
plot(u);
hold on;



u0=u;
u1=u;
u2=u;


for cout=0:dt:t
    
for i=4:N-1   
    u0(i)=(-2*u(i-3)+15*u(i-2)-60*u(i-1)+20*u(i)+30*u(i+1)-3*u(i+2))/(60*dx);
     %u0(i)=(-0.01027*u(i-3)+0.2180*u(i-2)-1.103*u(i-1)+0.6030*u(i)+0.282*u(i+1)+0.01027*u(i+2))/dx;
     %u0(i)=(-0.045*u(i-3)+0.308*u(i-2)-1.12*u(i-1)+0.45*u(i)+0.44*u(i+1)-0.038*u(i+2))/dx;
end
for i=4:N-1
    u1(i)=u(i)-dt*u0(i);
end

for i=4:N-1   
    u0(i)=(-2*u1(i-3)+15*u1(i-2)-60*u1(i-1)+20*u1(i)+30*u1(i+1)-3*u1(i+2))/(60*dx);
    %u0(i)=(-0.01027*u1(i-3)+0.2180*u1(i-2)-1.103*u1(i-1)+0.6030*u1(i)+0.282*u1(i+1)+0.01027*u1(i+2))/dx;
    %u0(i)=(-0.045*u1(i-3)+0.308*u1(i-2)-1.12*u1(i-1)+0.45*u1(i)+0.44*u1(i+1)-0.038*u1(i+2))/dx;
end
for i=4:N-1
    u2(i)=(3/4)*u(i)+(u1(i)-dt*u0(i))/4;
end

for i=4:N-1   
    u0(i)=(-2*u2(i-3)+15*u2(i-2)-60*u2(i-1)+20*u2(i)+30*u2(i+1)-3*u2(i+2))/(60*dx);
    %u0(i)=(-0.01027*u2(i-3)+0.2180*u2(i-2)-1.103*u2(i-1)+0.6030*u2(i)+0.282*u2(i+1)+0.01027*u2(i+2))/dx;
    %u0(i)=(-0.045*u2(i-3)+0.308*u2(i-2)-1.12*u2(i-1)+0.45*u2(i)+0.44*u2(i+1)-0.038*u2(i+2))/dx;
end
for i=4:N-1
    u(i)=u(i)/3+2*(u2(i)-dt*u0(i))/3;
end




u(1)=u(N-3);
u(2)=u(N-2);
u(3)=u(N-1);
u(N+1)=u(5);
u(N)=u(4);


end

figure;
plot(uu);
hold on;
plot(u,'--');
title('最高精度,t_total=20dt,局部放大');
   xlabel('x');
   ylabel('u');