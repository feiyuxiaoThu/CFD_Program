%求解线性波动方程
%Lax格式
Mx=100;
CFL=0.5;

a=1;
t_total=10; %0.1,1.0,10

dt_record=0.05;


Dx=0.5-(-0.5);
dx=Dx/Mx;

dt=dx*CFL/a;

Mt=ceil(t_total/dt);  %时间网格

Ux=zeros(1,Mx);
Uxx=zeros(1,Mx);
y=zeros(1,Mx);

x=zeros(1,Mx);

for i=1:1:Mx
    x(i)=-0.5+dx*(i-1);
      
end





%初始化
for i=1:1:Mx
    Ux(1,i)=fx(x(i));
    y(1,i)=fx(x(i)-rem(t_total,Dx));
end

Uxx(1,:)=Ux(1,:);

j=1;

for i=1:1:Mt
  
    u0=Uxx(1,Mx);
    u1=Uxx(1,1);
    Ux(1,1)=(0.5-CFL/2)*Uxx(1,2)+(0.5+CFL/2)*u0;
    Ux(1,Mx)=(0.5-CFL/2)*u1+(0.5+CFL/2)*Uxx(1,Mx-1);
    for ii=2:1:Mx-1
        
        Ux(1,ii)=(0.5-CFL/2)*Uxx(1,ii+1)+(0.5+CFL/2)*Uxx(1,ii-1);
    
    end
    
     
    
    
    Uxx(1,:)=Ux(1,:);
    
%     plot(Uxx(1,:));
%     hold on;
    
end

plot(Uxx);
  hold on;
 plot(y);
   %axis([-0.5 0.5 0 t_total 0 1]);
   title('Lax格式,t_total=10');
   xlabel('x');
   ylabel('t');


