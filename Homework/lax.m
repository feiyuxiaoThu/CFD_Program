%¼¤²¨¹Ü
%lax
clear;
dx=0.001;
dt=0.0001; %ÎÈ¶¨ÐÔ
CFL=0.1;
t_total=0.25;
L=1;
gamma=1.4;
x=-0.5:dx:0.5;
num=L/dx+1;
ro=zeros(1,num);
rro=zeros(1,num);
p=zeros(1,num); 
pp=zeros(1,num);
u=zeros(1,num);
uu=zeros(1,num);
e=zeros(1,num);
ee=zeros(1,num);

Lt=zeros(1,num);




for i=1:(num-1)/2 +1;
    u(i)=0.75;    
    ro(i)=1;
    p(i)=1;
end

for i=(num-1)/2+2:num
    u(i)=0;
    ro(i)=0.125;
    p(i)=0.1;
end

for i=1:num
    e(i)=p(i)/(ro(i)*(gamma-1))+0.5*u(i)*u(i);
    ee(i)=e(i);
    rro(i)=ro(i);
    pp(i)=p(i);
    uu(i)=u(i);
end

t=0;

while 1
    for i=1:num
        Lt(i)=abs(u(i))+(gamma*p(i)/ro(i))^0.5;
    end
    
   dt=CFL*dx/max(Lt);

    
    
    for i=2:num-1
        
         
             
        
        F1=0.5*(ro(i+1)*u(i+1)+ro(i)*u(i))-0.5*(dx/dt) *(ro(i+1)-ro(i));
        F2=0.5*(ro(i)*u(i)+ro(i-1)*u(i-1))-0.5*(dx/dt) *(ro(i)-ro(i-1));
        rro(i)=ro(i)-(dt/dx)*(F1-F2);
        
        FF1=0.5*(ro(i+1)*u(i+1)^2+p(i+1)+ro(i)*u(i)^2+p(i))-0.5* (dx/dt) *(ro(i+1)*u(i+1)-ro(i)*u(i));
        FF2=0.5*(ro(i)*u(i)^2+p(i)+ro(i-1)*u(i-1)^2+p(i-1))-0.5* (dx/dt) *(ro(i)*u(i)-ro(i-1)*u(i-1));
        uu(i)=(1/rro(i))*(ro(i)*u(i)-(dt/dx) *(FF1-FF2));
        
         FFF1=0.5*(u(i+1)*e(i+1)*ro(i+1)+p(i+1)*u(i+1)+u(i)*e(i)*ro(i)+p(i)*u(i))-0.5*(dx/dt) *(ro(i+1)*e(i+1)-ro(i)*e(i));
        FFF2=0.5*(u(i)*e(i)*ro(i)+p(i)*u(i)+u(i-1)*e(i-1)*ro(i-1)+p(i-1)*u(i-1))-0.5*(dx/dt) *(ro(i)*e(i)-ro(i-1)*e(i-1));
        ee(i)=(ro(i)*e(i)-(dt/dx) *(FFF1-FFF2))/rro(i);
        
        pp(i)=(ee(i)*rro(i)-0.5*rro(i)*uu(i)^2)*(gamma-1);
        
   
    end
    
%         F1=0.5*(ro(2)*u(2)+ro(1)*u(1))-0.5*dx/dt *(ro(2)-ro(1));
%         F2=0.5*(ro(1)*u(1)+ro(num)*u(num))-0.5*dx/dt *(ro(1)-ro(num));
%         rro(1)=ro(1)-dt/dx*(F1-F2);
%         FF1=0.5*(ro(2)*u(2)^2+Cpro(2)*ro(2)^gamma+ro(1)*u(1)^2+Cpro(1)*ro(1)^gamma)-0.5*dx/dt *(ro(2)*u(2)-ro(1)*u(1));
%         FF2=0.5*(ro(1)*u(1)^2+Cpro(1)*ro(1)^gamma+ro(num)*u(num)^2+Cpro(num)*ro(num)^gamma)-0.5*dx/dt *(ro(1)*u(1)-ro(num)*u(num));
%         uu(1)=(1/ro(1))*(ro(1)*u(1)-dt/dx *(FF1-FF2));
%         
%         F1=0.5*(ro(1)*u(1)+ro(num)*u(num))-0.5*dx/dt *(ro(1)-ro(num));
%         F2=0.5*(ro(num)*u(num)+ro(num-1)*u(num-1))-0.5*dx/dt *(ro(num)-ro(num-1));
%         rro(num)=ro(num)-dt/dx*(F1-F2);
%         FF1=0.5*(ro(1)*u(1)^2+Cpro(1)*ro(1)^gamma+ro(num)*u(num)^2+Cpro(num)*ro(num)^gamma)-0.5*dx/dt *(ro(1)*u(1)-ro(num)*u(num));
%         FF2=0.5*(ro(num)*u(num)^2+Cpro(num)*ro(num)^gamma+ro(num-1)*u(num-1)^2+Cpro(num-1)*ro(num-1)^gamma)-0.5*dx/dt *(ro(num)*u(num)-ro(num-1)*u(num-1));
%         uu(i)=(1/ro(num))*(ro(num)*u(num)-dt/dx *(FF1-FF2));
    
  
    for i=2:num-1
        ro(i)=rro(i);
        u(i)=uu(i);
        e(i)=ee(i);
        p(i)=pp(i);
    end
   
    t=t+dt;
    if t>=t_total
        break;
    end
end


figure(1)
plot(x,u,'-b','LineWidth',2)
xlabel('Position')
ylabel('Velocity')

figure(2)
plot(x,p,'-b','LineWidth',2)
xlabel('Position')
ylabel('Pressure')

figure(3)
plot(x,ro,'-b','LineWidth',2)
xlabel('Position')
ylabel('Density') 
    