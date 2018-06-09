%求解线性波动方程
%Warming-Beam格式
clear all;
Mx=256;
N=Mx;
k0=5;
CFL=0.1;
cc=CFL*CFL;
a=1;

dt_record=0.005;

Dx=0.5-(-0.5);
dx=Dx/Mx;

dt=dx*CFL/a;
t_total=20*dt;

Mt=ceil(t_total/dt);  %时间网格

Ux=zeros(1,Mx+1);
Uxx=zeros(1,Mx+1);
y=zeros(1,Mx+1);

x=zeros(1,Mx+1);

for i=1:1:Mx+1
    x(i)=-0.5+dx*(i-1);
end


Ux(:)=1;



for i=1:N+1
    for j=1:N/4
        E=(j/k0)^4*exp(-2*(j/k0)^2);
    Ux(i)=Ux(i)+0.1*(E^0.5*sin(2*pi*j*(x(i)+rand(1,1))));
    end
end

% for i=1:N+1
%     Ux(i)=3*sin(2*pi*x(i))+4*sin(100*pi*x(i))+5*sin(50*pi*x(i));
% end
% for i=1:N+1
%     UUx(i)=3*sin(2*pi*(x(i)-t_total))+4*sin(100*pi*(x(i)-t_total))+5*sin(50*pi*(x(i)-t_total));
% end
plot(Ux);
hold on;

Uxx(1,:)=Ux(1,:);

for i=1:1:Mt
    
    u0=Uxx(1,Mx-2);
    um1=Uxx(1,Mx-1);
    Ux(1,1)=(1-1.5*CFL+cc/2)*Uxx(1,1)+(2*CFL-cc)*u0+(cc/2-CFL/2)*um1;
    Ux(1,2)=(1-1.5*CFL+cc/2)*Uxx(1,2)+(2*CFL-cc)*Uxx(1,1)+(cc/2-CFL/2)*u0;
    for ii=3:1:Mx;
        
        Ux(1,ii)=(1-1.5*CFL+cc/2)*Uxx(1,ii)+(2*CFL-cc)*Uxx(1,ii-1)+(cc/2-CFL/2)*Uxx(1,ii-2);
    
    end
    
    
    
    
    Uxx(1,:)=Ux(1,:);
    
%     plot(Uxx(1,:));
%     hold on;
    
end

    
 
 plot(Uxx,'--');
     
   %axis([-0.5 0.5 0 t_total 0 1]);
   title('Warming-Beam格式,t_total=20dt');
   xlabel('x');
   ylabel('u');


