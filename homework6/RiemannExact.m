function [] = RiemannExact(p1,rho1,u1,p4,rho4,u4,tol)

%This code gives the exact solution of Riemann's 1-D problem in the context
%of the unsteady 1-D shock tube,by means of the Newton-Raphson's method. 
%After having set the inputs to the left and right
%gases' variables, one can easily check the waves generated from the
%discontinuity; the program returns in output the type of the solution,
%which can be for example an RCS (meaning that an expansion has generated
%on the first family and a shock on the second family). The user then has
%the chance to get the most important variables (sound speed, pressure,
%velocity) plotted versus the shock tube's length coordinate. NB: all
%interactions are neglected. 
%
%INPUTS:
%
%p1: pressure of the gas in the left part of the shock tube
%rho1: density of the gas in the left part of the shock tube
%u1: velocity of the particles in the left part of the shock tube
%p4: pressure of the gas in the right part of the shock tube
%rho4: density of the gas in the right part of the shock tube
%u4: velocity of the particles in the right part of the shock tube
%tol: tolerance of solution
%
%Virginia Notaro


iter=0;
gamma=1.4;
delta=(gamma-1)/2;
alfa=2*gamma/(gamma-1);

a1=sqrt(gamma*p1/rho1);
a4=sqrt(gamma*p4/rho4);

zeta=(p1/p4)^(1/alfa)*(a1/a4); %Entropy difference between gas 1 and gas 4

if u4-u1>a4/delta+a1/delta
    disp('Warning: Solution does not exist!')
else
    
    %First approximation velocity
    ustar=(zeta*(a1+delta*u1)-a4+delta*u4)/((1+zeta)*delta);
    u=ustar;
    p2=1;
    p3=0.5;
    
    while abs(p2-p3)>tol
        
    %Checking wave on first family
    if u>u1  %R
    a2=a1+delta*(u1-u);
    p2=p1*(a2/a1)^alfa;
    dp2=-gamma*p2/a2;
    
    elseif u<u1 %S
    V1=(gamma+1)/4*(u1-u)+sqrt(((gamma+1)/4)^2*(u1-u)^2+a1^2);
    M1=V1/a1;
    p2=p1*(1+(2*gamma)/(gamma+1)*(M1^2-1));
    a2=a1*sqrt((gamma+1+(gamma-1)*(p2/p1))/(gamma+1+(gamma-1)*(p1/p2)));
    ws=u1-V1;
    M2=(u-ws)/a2;
    dp2=-2*gamma*(p1/a1)*((abs(M1))^3/(1+M1^2));
    
    
    
    end

    %Checking wave on the second family
    if u>u4 %S
    V4=(gamma+1)/4*(u-u4)+sqrt(((gamma+1)/4)^2*(u-u4)^2+a4^2);
    M4=V4/a4;
    p3=p4*(1+(2*gamma)/(gamma+1)*(M4^2-1));
    V3=V4*(1+delta*M4^2)/((gamma+1)/2*M4^2);
    a3=a4*sqrt((gamma+1+(gamma-1)*(p3/p4))/(gamma+1+(gamma-1)*(p4/p3)));
    wd=u4-V4;
    M3=(u-wd)/a3;
    dp3=2*gamma*(p4/a4)*((abs(M4))^3/(1+M4^2));
    
    elseif u<u4 %R
    a3=a4+delta*(u-u4);
    p3=p4*(a3/a4)^alfa;
    dp3=gamma*p3/a3;
    
    end
    
    %New velocity and update of the iteration
    u=u-(p2-p3)/(dp2-dp3);
    iter=iter+1;
    
    %Checking number of iterations
    if iter>45000
        disp('The problem can not be solved using this method!')
        break
    end
    end
    
    %End of cycle and density on the right (3) and left (2) side of the
    %contact discontinuity
    
    rho2=gamma*p2/a2^2;
    rho3=gamma*p3/a3^2;
 
    
    %Type of the solution
  
    
    if ((u>u1)&&(abs(p2-p1)>tol)&&(abs(p2-p4)>tol))
        if (u<u4)
            tipo=('RCR');
            disp('The solution is RCR')
            
        else
            tipo=('RCS');
            disp('The solution is RCS')
            disp('wd= '),disp(wd)
        end
    elseif ((u<u1)&&(abs(p2-p1)>tol)&&(abs(p2-p4)>tol))
        if (u<u4)
            tipo=('SCR');
            disp('The solution is SCR')
            disp('ws= '),disp(ws)
        else
            tipo=('SCS');
            disp('The solution is SCS')
            disp('ws= '),disp(ws)
            disp('wd= '),disp(wd)
        end
    elseif ((abs(u-u1)<tol)&&(abs(p2-p1)<tol)&&(abs(u-u4)>tol))
        if (u<u4)
            tipo=('NCR');
            disp('The solution is NCR')
            u=u1;
            p2=p1; 
            p3=p1;
            a2=a1; 
           
            rho2=rho1; 
            
                      
        else
            tipo=('NCS');
            disp('The solution is NCS')
            u = u1;
            p2 = p1; 
            p3 = p1;
            a2 = a1; 
            disp('wd= '),disp(wd)
            
            rho2 = rho1; 
            
           
        end
    elseif ((abs(u-u1)>tol)&&(abs(p2-p4)<tol)&&(abs(u-u4)<tol))
        if (u>u1)
            tipo=('RCN');
            disp('The solution is RCN')
            u = u4;
            p2 = p4; 
            p3 = p4;
            a3 = a4;
            rho3 = rho4;
           
        else
            tipo=('SCN');
            disp('The solution is SCN')
            u = u4;
            p2 = p4; 
            p3 = p4;
            disp('ws= '),disp(ws)
             
            a3 = a4;
            rho3 = rho4;
            
        end
    else
        tipo=('NCN');
        disp('There isn''t any wave! NCN')
        u = u1;
        p2 = p1; p3 = p1;
        a2 = a1; a3 = a1;
        rho2 = rho1; rho3 = rho1;
              
    end
    
    disp('u= '),disp(u)
    disp('p2= '),disp(p2)
    disp('p3= '),disp(p3)
    disp('rho2= '),disp(rho2)
    disp('rho3= '),disp(rho3)
    
    
end
    
    disp('The number of iterations is: '),disp(iter)
    
    
%Graphical representation of pressure, velocity and sound speed versus the
%position x in the shock tube. 

t=input('Insert the time which you want to calculate the solution at (check the velocities'' order of magnitude first!), t= ');
x0=input('Insert the lower bound of the visualization domain, x0= ');
xf=input('Insert the upper bound of the visualization domain, xf= ');
npoints=input('insert the number of points which you want to calculate the solution in, npoints= ');

x=linspace(x0,xf,npoints);

A = zeros(npoints,1);   
P = zeros(npoints,1); 
U = zeros(npoints,1);
rofinal=zeros(npoints,1);


x1=(u1-a1)*t;
x2=(u-a2)*t;
x3=u*t;
x4=(u+a3)*t;
x5=(u4+a4)*t;

%Analytical solution by means of the characteristic equations

if tipo==('RCR')
    
    for k=1:npoints
        if x(k)<=x1
            
            U(k)=u1;
            P(k)=p1;
            A(k)=a1;
            
        elseif x(k)<=x2 && x(k)>=x1
            
            U(k)=1/(1+delta)*(a1+delta*u1+(x(k))/t);
            A(k)=1/(1+delta)*(a1)+delta/(delta+1)*(u1)-delta/(1+delta)*(x(k))/t;
            P(k)=p1.*(A(k)/a1).^(alfa);
            
        elseif x(k)>=x2 && x(k)<=x3
            
            U(k)=u;
            A(k)=a2;
            P(k)=p2;
            
        elseif x(k)>=x3 && x(k)<=x4
            
            U(k)=u;
            A(k)=a3;
            P(k)=p3;
            
        elseif x(k)>=x4 && x(k)<=x5
            
            U(k)=1/(1+delta)*((x(k))/t-a4+delta*u4);
            A(k)=1/(1+delta)*(a4)-delta/(1+delta)*(u4)+delta/(1+delta)*(x(k))/t;
            P(k)=p4.*(A(k)./a4).^(alfa);
            
        elseif x(k)>=x5
            
            U(k)=u4;
            P(k)=p4;
            A(k)=a4;
        end
    end
elseif tipo==('RCS')
    
    xUD=abs(wd)*t;
    
    for k=1:npoints
        if x(k)<=x1
            
            U(k)=u1;
            P(k)=p1;
            A(k)=a1;
            
        elseif x(k)<=x2 && x(k)>=x1
            
            U(k)=1/(1+delta)*(a1+delta*u1+(x(k))/t);
            A(k)=1/(1+delta)*(a1)+delta/(delta+1)*(u1)-delta/(1+delta)*(x(k))/t;
            P(k)=p1.*(A(k)./a1).^(alfa);
            
        elseif x(k)>=x2 && x(k)<=x3
            
            U(k)=u;
            A(k)=a2;
            P(k)=p2;
            
        elseif x(k)>=x3 && x(k)<=xUD
            
            U(k)=u;
            P(k)=p3;
            A(k)=a3;
            
        elseif x(k)>=xUD
            
            U(k)=u4;
            P(k)=p4;
            A(k)=a4;
        end
    end
    
elseif tipo==('SCR')
    
    xUS=-abs(ws)*t;
    
    for k=1:npoints;
        
        if x(k)<=xUS
            U(k)=u1;
            P(k)=p1;
            A(k)=a1;
            
        elseif x(k)>=xUS && x(k)<=x3
            
            U(k)=u;
            P(k)=p2;
            A(k)=a2;
            
        elseif x(k)>=x3 && x(k)<=x4
            
            U(k)=u;
            A(k)=a3;
            P(k)=p3;
            
        elseif x(k)>=x4 && x(k)<=x5
            
            U(k)=1/(1+delta)*((x(k))/t-a4+delta*u4);
            A(k)=1/(1+delta)*(a4)-delta/(delta+1)*(u4)+delta/(1+delta)*(x(k))/t;
            P(k)=p4.*(A(k)./a4).^(alfa);
            
        elseif x(k)>=x5
            
            U(k)=u4;
            P(k)=p4;
            A(k)=a4;
        end
    end
    
elseif tipo==('SCS')
    
    xUS=-abs(ws)*t;
    
    xUD=abs(wd)*t;
    
    for k=1:npoints
        if x(k)<=xUS
            U(k)=u1;
            P(k)=p1;
            A(k)=a1;
            
        elseif x(k)>=xUS && x(k)<=x3
            
            U(k)=u;
            P(k)=p2;
            A(k)=a2;
            
        elseif x(k)>=x3 && x(k)<=xUD
            
            U(k)=u;
            P(k)=p3;
            A(k)=a3;
            
        elseif x(k)>=xUD
            
            U(k)=u4;
            P(k)=p4;
            A(k)=a4;
        end
    end
    
elseif tipo==('NCR')
    
    for k=1:npoints
        if x(k)<=x3
            
            U(k)=u1;
            A(k)=a1;
            P(k)=p1;
            
        elseif x(k)>=x3 && x(k)<=x4
            
            U(k)=u;
            A(k)=a3;
            P(k)=p3;
            
        elseif x(k)>=x4 && x(k)<=x5
            
            U(k)=1/(1+delta)*((x(k))/t-a4+delta*u4);
            A(k)=1/(1+delta)*(a4)-delta/(delta+1)*(u4)+delta/(1+delta)*(x(k))/t;
            P(k)=p4.*(A(k)./a4).^(alfa);
            
        elseif x(k)>=x5
            
            U(k)=u4;
            P(k)=p4;
            A(k)=a4;
            
        end
    end
    
elseif tipo==('NCS')
    
    xUD=abs(wd)*t;
    
    for k=1:npoints
        if x(k)<=x3
            
            U(k)=u1;
            A(k)=a1;
            P(k)=p1;
            
        elseif x(k)>=x3 && x(k)<=xUD
            
            U(k)=u;
            P(k)=p3;
            A(k)=a3;
            
        elseif x(k)>=xUD
            
            U(k)=u4;
            P(k)=p4;
            A(k)=a4;
            
        end
    end
    
elseif tipo==('RCN')
    
    for k=1:npoints
        if x(k)<=x1
            
            U(k)=u1;
            P(k)=p1;
            A(k)=a1;
            
        elseif x(k)>=x1 && x(k)<=x2
            
            U(k)=1/(1+delta)*(a1+delta*u1+(x(k))/t);
            A(k)=1/(1+delta)*(a1)+delta/(delta+1)*(u1)-delta/(1+delta)*(x(k))/t;
            P(k)=p1.*(A(k)./a1).^(alfa);
            
        elseif x(k)>=x2 && x(k)<=x3
            
            U(k)=u;
            A(k)=a2;
            P(k)=p2;
            
        elseif x(k)>=x3
            
            U(k)=u4;
            A(k)=a4;
            P(k)=p4;
            
        end
    end
    
elseif tipo==('SCN')
    
    xUS=-abs(ws)*t;
    
    for k=1:npoints
        
        if x(k)<=xUS
            
            U(k)=u1;
            P(k)=p1;
            A(k)=a1;
            
        elseif x(k)>=xUS && x(k)<=x3
            
            U(k)=u;
            P(k)=p2;
            A(k)=a2;
            
        elseif x(k)>=x3
            
            U(k)=u4;
            A(k)=a4;
            P(k)=p4;
            
        end
    end
    
elseif tipo==('NCN')
    
    for k=1:npoints
        
        U(k)=u1;
        P(k)=p1;
        A(k)=a1;
        
    end
end

for i=1:npoints
rofinal(i)=A(i)^2*1.4*P(i);
end

%Data plot
    
grid on

figure(1)
plot(x,U,'-b','LineWidth',2)
xlabel('Position')
ylabel('Velocity')

figure(2)
plot(x,P,'-b','LineWidth',2)
xlabel('Position')
ylabel('Pressure')

figure(3)
plot(x,rofinal,'-b','LineWidth',2)
xlabel('Position')
ylabel('Density') 

end

