% Jamson格式
%激波管
%lax
clear;
dx=0.001;
dt=0.01; %稳定性
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

%系数
v=zeros(1,4);
vv=zeros(1,4);
e2=zeros(1,2);
e4=zeros(1,2);

k2=3; %0.5-0.6
k4=0.01;

EE=zeros(1,2);
EEE=zeros(1,2);


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
    
    
    for i=4:num-3
        
         L0=abs(u(i-1))+(gamma*p(i-1)/ro(i-1))^0.5;
          L1=abs(u(i))+(gamma*p(i)/ro(i))^0.5;
           L2=abs(u(i+1))+(gamma*p(i+1)/ro(i+1))^0.5;
           
           if (L0>L1)
               LL1=L0;
           else
               LL1=L1;
           end
           
           if (L1>L2)
               LL2=L1;
           else
               LL2=L2;
           end
           
           for j=0:3
               v(j+1)=abs(p(i+3-j)-2*p(i+2-j)+p(i+1-j))/abs(p(i+3-j)+2*p(i+2-j)+p(i+1-j));
           end
           for j=0:3
               vv(j+1)=abs(p(i+2-j)-2*p(i+1-j)+p(i-j))/abs(p(i+2-j)+2*p(i+1-j)+p(i-j));
           end
           
           e2(1)=max(v)*k2;
           e2(2)=max(vv)*k2;
           
           EE(1)=0;
           EE(2)=k4-e2(1);
           EEE(1)=0;
           EEE(2)=k4-e2(2);
           
           e4(1)=max(EE);
           e4(2)=max(EEE);
           
           
           
           
             
        
        F1=0.5*(ro(i+1)*u(i+1)+ro(i)*u(i))-LL1*e2(1) *(ro(i+1)-ro(i))+LL1*e4(1)*(ro(i+2)-3*ro(i+1)+3*ro(i)-ro(i-1));
        F2=0.5*(ro(i)*u(i)+ro(i-1)*u(i-1))-LL2*e2(2) *(ro(i)-ro(i-1))+LL2*e4(2)*(ro(i+1)-3*ro(i)+3*ro(i-1)-ro(i-2));
        rro(i)=ro(i)-(dt/dx)*(F1-F2);
        
        FF1=0.5*(ro(i+1)*u(i+1)^2+p(i+1)+ro(i)*u(i)^2+p(i))-LL1*e2(1) *(ro(i+1)*u(i+1)-ro(i)*u(i))+LL1*e4(1)*(ro(i+2)*u(i+2)-3*ro(i+1)*u(i+1)+3*ro(i)*u(i)-ro(i-1)*u(i-1));
        FF2=0.5*(ro(i)*u(i)^2+p(i)+ro(i-1)*u(i-1)^2+p(i-1))-LL2*e2(2) *(ro(i)*u(i)-ro(i-1)*u(i-1))+LL2*e4(2)*(ro(i+1)*u(i+1)-3*ro(i)*u(i)+3*ro(i-1)*u(i-1)-ro(i-2)*u(i-2));
        uu(i)=(1/rro(i))*(ro(i)*u(i)-(dt/dx) *(FF1-FF2));
        
         FFF1=0.5*(ro(i+1)*u(i+1)*e(i+1)+p(i+1)*u(i+1)+ro(i)*u(i)*e(i)+p(i)*u(i))-LL1*e2(1) *(ro(i+1)*e(i+1)-ro(i)*e(i))+LL1*e4(1)*(ro(i+3)*e(i+2)-3*ro(i+1)*e(i+1)+3*ro(i)*e(i)-ro(i-1)*e(i-1));
        FFF2=0.5*(u(i)*e(i)*ro(i)+p(i)*u(i)+u(i-1)*e(i-1)*ro(i-1)+p(i-1)*u(i-1))-LL2*e2(2) *(ro(i)*e(i)-ro(i-1)*e(i-1))+LL2*e4(2)*(ro(i+1)*e(i+1)-3*ro(i)*e(i)+3*ro(i-1)*e(i-1)-ro(i-2)*e(i-2));
        ee(i)=(ro(i)*e(i)-(dt/dx) *(FFF1-FFF2))/rro(i);
        
        pp(i)=(rro(i)*ee(i)-0.5*rro(i)*uu(i)^2)*(gamma-1);
        
   
    end
    
    
  
    for i=4:num-3
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

    