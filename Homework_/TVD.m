%TVD
clear;
t=10;%1.0,10.0
Mx=100;
L=0.5;
dx=2*L/Mx;
numx=Mx+1;
x=zeros(1,numx);
for i=1:numx
    x(i)=-0.5+dx*(i-1);
end
CFL=0.5;
a=1;


u0=zeros(1,numx+3);
u=zeros(1,numx);

for i=26:76
    u(i)=1;
end

U=u;

for i=1:numx
    u0(i+2)=u(i);
end

u0(1)=u(numx-1);
u0(2)=u(numx);
u0(numx+3)=u(1);
 %plot(u);
 
 uu=zeros(1,numx);
 
 dt=CFL*dx/a;
 
 for t_total=0:dt:t
     for j=3:numx+2
         
         D1=0;
         D0=0;
         if (u0(j+1)-u0(j))*(u0(j)-u0(j-1))>0
             if abs(u0(j+1)-u0(j))<abs(u0(j)-u0(j-1))
             D1=sign(u0(j+1)-u0(j))*abs(u0(j+1)-u0(j));
             else
             D1=sign(u0(j+1)-u0(j))*abs(u0(j)-u0(j-1));
             end
         end
         
         if (u0(j)-u0(j-1))*(u0(j-1)-u0(j-2))>0
             if abs(u0(j)-u0(j-1))<abs(u0(j-1)-u0(j-2))
             D0=sign(u0(j)-u0(j-1))*abs(u0(j)-u0(j-1));
             else
             D0=sign(u0(j)-u0(j-1))*abs(u0(j-1)-u0(j-2));
             end
         end
           D0=D0/dx;
           D1=D1/dx;
           
           f1=a*(u0(j)+0.5*D1*(dx-a*dt));
           f0=a*(u0(j-1)+0.5*D0*(dx-a*dt));
         
        
         
         uu(j-2)=u0(j)+(dt/dx)*(f0-f1);
         
     end
     
     for i=1:numx
        u0(i+2)=uu(i);
     end

    u0(1)=uu(numx-1);
    u0(2)=uu(numx);
    u0(numx+3)=uu(1);
     
     
 end
 
 for i=1:numx
     u(i)=u0(i+2);
 end
 
 
 
 plot(x,u);
 hold on;
 plot(x,U);
 xlabel('x');
 ylabel('u');
 
 