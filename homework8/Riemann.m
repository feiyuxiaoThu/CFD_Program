%Riemann Problem
clear;

L=1; %求解区域
node=50; %网格数
dx=2*L/node;
dt=0.0001;
t_tol=0.01;
x=-L:dx:(L-dx);

CFL=0.1;

ga=1.4; %\gamma

ro=zeros(1,node);
u=zeros(1,node);
e=zeros(1,node);

p=zeros(1,node);


%初值
for i=1:25
    ro(i)=1;
    u(i)=0;
    p(i)=1;
    e(i)=p(i)/((ga-1)*ro(i));
end

for i=26:50
    ro(i)=0.125;
    u(i)=0;
    p(i)=0.1;
    e(i)=p(i)/((ga-1)*ro(i));
end

rro=ro;
uu=u;
ee=e;
pp=p;


for t=0:dt:t_tol

 A=zeros(3,3,node-3);
%零阶重构

%线性化
for j=1:node-3
    i=j+1;
    U=zeros(1,3);
    U(1)=ro(i);
    U(2)=ro(i)*u(i);
    U(3)=ro(i)*e(i);
    
   
    A(1,2,j)=1;
    A(2,1,j)=(ga-3)*u(i)*u(i)/2;
    A(2,2,j)=(3-ga)*u(i);
    A(2,3,j)=ga-1;
    A(3,1,j)=(ga-1)*u(i)^3-ga*u(i)*e(i);
    A(3,2,j)=-1.5*(ga-1)*u(i)*u(i)+ga*e(i);
    A(3,3,j)=ga*u(i);
end
 L=zeros(3,3,node-3);
 R=zeros(3,3,node-3);
 I=zeros(3,3,node-3);
    
for j=1:node-3
    
    [V,D]=eig(A(:,:,j));
   
    L(:,:,j)=inv(V);
    R(:,:,j)=V;
    I(:,:,j)=D;
      
end

F=zeros(3,node-3);

Lt=zeros(1,node-3);

for i=1:node-3
    Lt(i)=abs(u(i))+(ga*p(i)/ro(i))^0.5;
end

dt=CFL*dx/max(Lt);

%重构
for j=1:node-3
    U1=zeros(1,3);
    U2=zeros(1,3);
    U3=zeros(1,3);
    U4=zeros(1,3);
    
    i=j+1;
    
    U1(1)=ro(i-1);
    U1(2)=ro(i-1)*u(i-1);
    U1(3)=ro(i-1)*e(i-1);
    
    U2(1)=ro(i);
    U2(2)=ro(i)*u(i);
    U2(3)=ro(i)*e(i);
    
    U3(1)=ro(i+1);
    U3(2)=ro(i+1)*u(i+1);
    U3(3)=ro(i+1)*e(i+1);
    
    U4(1)=ro(i+2);
    U4(2)=ro(i+2)*u(i+2);
    U4(3)=ro(i+2)*e(i+2);
    
    a1=L(1,:,j)*(U3-U2)';
    b1=L(1,:,j)*(U2-U1)';
    if a1*b1>0
        if abs(a1)>abs(b1)
            D1=sign(a1)*abs(b1)/dx;
        else
            D1=sign(a1)*abs(a1)/dx;
        end
    else
        D1=0;
    end
    a2=L(2,:,j)*(U3-U2)';
    b2=L(2,:,j)*(U2-U1)';
    if a2*b2>0
        if abs(a2)>abs(b2)
            D2=sign(a2)*abs(b2)/dx;
        else
            D2=sign(a2)*abs(a2)/dx;
        end
    else
        D2=0;
    end
    a3=L(3,:,j)*(U3-U2)';
    b3=L(3,:,j)*(U2-U1)';
    if a3*b3>0
        if abs(a3)>abs(b3)
            D3=sign(a3)*abs(b3)/dx;
        else
            D3=sign(a3)*abs(a3)/dx;
        end
    else
        D3=0;
    end
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aa1=L(1,:,j)*(U4-U3)';
    bb1=L(1,:,j)*(U3-U2)';
    if aa1*bb1>0
        if abs(aa1)>abs(bb1)
            DD1=sign(aa1)*abs(bb1);
        else
            DD1=sign(aa1)*abs(aa1);
        end
    else
        DD1=0;
    end
     aa2=L(1,:,j)*(U4-U3)';
    bb2=L(1,:,j)*(U3-U2)';
    if aa2*bb2>0
        if abs(aa2)>abs(bb2)
            DD2=sign(aa1)*abs(bb2);
        else
            DD2=sign(aa2)*abs(aa2);
        end
    else
        DD2=0;
    end
     aa3=L(1,:,j)*(U4-U3)';
    bb3=L(1,:,j)*(U3-U2)';
    if aa3*bb3>0
        if abs(aa3)>abs(bb3)
            DD3=sign(aa3)*abs(bb3);
        else
            DD3=sign(aa3)*abs(aa3);
        end
    else
        DD3=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WL=zeros(1,3);
    WR=zeros(1,3);
    
    WL(1)=L(1,:,j)*U2'+0.5*D1*(dx-I(1,1,j)*dt);
    WL(2)=L(2,:,j)*U2'+0.5*D2*(dx-I(2,2,j)*dt);
    WL(3)=L(3,:,j)*U2'+0.5*D3*(dx-I(3,3,j)*dt);
    
    WR(1)=L(1,:,j)*U3'+0.5*DD1*(dx-I(1,1,j)*dt);
    WR(2)=L(2,:,j)*U3'+0.5*DD2*(dx-I(2,2,j)*dt);
    WR(3)=L(3,:,j)*U3'+0.5*DD3*(dx-I(3,3,j)*dt);
    
%     WW=zeros(1,3);
%     WW(1)=0.5*(1+sign(I(1,1,j)))*WL(1)+0.5*(1-sign(I(1,1,j)))*WR(1);
%     WW(2)=0.5*(1+sign(I(2,2,j)))*WL(2)+0.5*(1-sign(I(2,2,j)))*WR(2);
%     WW(3)=0.5*(1+sign(I(3,3,j)))*WL(3)+0.5*(1-sign(I(3,3,j)))*WR(3);
%    UL=zeros(1,3);
%    UR=zeros(1,3);
   UL=R(:,:,j)*WL';
   UR=R(:,:,j)*WR';
   
   UL(2)=UL(2)/UL(1);
   if UL(2)<0.000001
       UL(2)=0;
   end
   UL(3)=UL(3)/UL(1);
   
   UR(2)=UR(2)/UR(1);
   if UR(2)<0.000001
       UR(2)=0;
   end
   UR(3)=UR(3)/UR(1);
   
   ro_=(UL(1)*UR(1))^0.5;
   u_=(UL(1)^0.5*UL(2)+UR(1)^0.5*UR(2))/(UL(1)^0.5+UR(1)^0.5);
   HL=UL(3)+0.5*UL(2)*UL(2);
   HR=UR(3)+0.5*UR(2)*UR(2);
   H_=(UL(1)^0.5*HL+UR(1)^0.5*HR)/(UL(1)^0.5+UR(1)^0.5);
   
   c_2=(ga-1)*(H_-0.5*u_*u_);
   p_=c_2*ro_/ga;
   
   e_=H_-0.5*u_*u_;
   
   F(1,j)=ro_*u_;
   F(2,j)=ro_*u_*u_+p_;
   F(3,j)=(ro_*e_+p_)*u_;
    
    
    
end

%时间相关离散
%Runge-Kutta法
for i=1:node-4
    j=i+1;
    
    rro(j)=dt*(F(1,i)-F(1,i+1))/dx+ro(j);
    uu(j)=dt*(F(2,i)-F(2,i+1))/dx+u(j);
    ee(j)=dt*(F(3,i)-F(3,i+1))/dx+e(j);
    
    pp(j)=(ga-1)*ee(j)*rro(j);
    
end


ro=rro;
p=pp;
e=ee;
u=uu;


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
    