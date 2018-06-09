%Simulating a Lid-Driven Cavity flow with a lid moving at speed U
...Flow considered is laminar and of an incompressible Newtonian fluid
...Pressure correction method is used on a staggered grid
...to solve the Navier Stokes equations

%%
clear;
%Specifying Parameters
nx=20;                        %Number of steps in space(x)
ny=20;                        %Number of steps in space(y)
T=5;                          %Final time
dt=1e-2;                      %Width of time step
nt=ceil(T/dt);                %Number of time steps
dt=T/nt;                      %Corrected time step
L=1;                          %Length of cavity
H=1;                          %Height of cavity
Re=1e3;                       %Flow Reynolds number
dx=L/nx;                  %Width of space step(x)
dy=H/ny;                  %Width of space step(y)
x=0:dx:L;                     %Range of x(0,L) and specifying grid points
y=0:dy:H;                     %Range of y(0,H) and specifying grid points
u=zeros(nx+3,ny+3);             %Preallocating u
v=zeros(nx+3,ny+3);             %Preallocating v
p=zeros(nx+3,ny+3);               %Preallocating p
pex=zeros(nx+3,ny+3);          %ÐÞÕýÁ¿
ppex=pex;


U=(1e-6)*Re;

UN=1;   VN=0;      
US=0;   VS=0;
UE=0;   VW=0;                  
UW=0;   VE=0;

%%
%Boundary conditions
u(1,:)=2*UW-u(2,:); u(end,:)=2*UE-u(end-1,:); u(:,1)=2*US-u(:,2); u(:,end)=2*UN-u(:,end-1);
v(1,:)=2*VW-v(2,:); v(end,:)=2*VE-v(end-1,:); v(:,1)=2*VS-v(:,2); v(:,end)=2*VN-v(:,end-1);

%%
%%
%Evaluating pressure and velocity field at each time step
uu=u;
vv=v;
nt=5;
for it=0:dt:1
   % velocity of x
  
   for i=2:nx+2
       for j=2:ny+2
          ac=dy*(0.25*(u(i+1,j)-u(i-1,j))+2/(Re*dx))+dx*(0.25*(v(i,j)+(v(i+1,j)-v(i,j-1)-v(i+1,j-1))+2/(Re*dy)))+dx*dy/dt;
          aright=dx*(0.25*(v(i,j)+v(i+1,j))-1/(Re*dy));
          aleft=-dx*(0.25*(v(i,j-1)+v(i+1,j-1))+1/(Re*dy));
          adown=-dy*(0.25*(u(i,j)+u(i-1,j))+1/(Re*dx));
          aup=dy*(0.25*(u(i+1,j)+u(i,j))-1/(Re*dx));   
          
          Ab=-u(i,j)*dx*dy/dt;
          
          uu(i,j)=((p(i,j)-p(i+1,j))*dy-Ab-aleft*uu(i,j-1)-adown*uu(i-1,j)-aright*u(i,j+1)-aup*u(i+1,j))/ac;
       end
   end
  
   
   
   % velocity of y
   
   for i=2:nx+2
       for j=2:ny+2
           bc=dx*(0.25*(v(i,j+1)-v(i,j-1))+2/(Re*dy))+dy*(0.25*(u(i,j+1)+u(i,j)-u(i-1,j+1)-u(i-1,j))+2/(Re*dx));
           bright=dx*(0.25*(v(i,j)+v(i,j+1))-1/(Re*dy));
           bleft=-dx*(0.25*(v(i,j-1)+v(i,j))+1/(Re*dy));
           bdown=-dy*(0.25*(u(i-1,j+1)+u(i-1,j))+1/(Re*dx));
           bup=dy*(0.25*(u(i,j+1)+u(i,j))-1/(Re*dx));    
           
           Bb=-v(i,j)*dx*dy/dt;
           
           vv(i,j)=((p(i,j)-p(i,j+1))*dy-Bb-bleft*vv(i,j-1)-bdown*vv(i-1,j)-bright*v(i,j+1)-bup*v(i+1,j))/bc;
       end
   end
   
    
   
   
   %Boundary conditions
    u(1,:)=2*UW-u(2,:); u(end,:)=2*UE-u(end-1,:); u(:,1)=2*US-u(:,2); u(:,end)=2*UN-u(:,end-1);
    v(1,:)=2*VW-v(2,:); v(end,:)=2*VE-v(end-1,:); v(:,1)=2*VS-v(:,2); v(:,end)=2*VN-v(:,end-1);
    
    
    % p
    %boundary conditions for p
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p(1,:)=p(2,:);p(end,:)=p(end-1,:);
    p(:,1)=p(:,2);p(:,end)=p(:,end-1);
    pp=p;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ex=zeros(nx+3,ny+3);
    
    val=1;
    
    while val>0.001
        for i=2:nx+2
            for j=2:ny+2
                ac=dy*(0.25*(u(i+1,j)-u(i-1,j))+2/(Re*dx))+dx*(0.25*(v(i,j)+(v(i+1,j)-v(i,j-1)-v(i+1,j-1))+2/(Re*dy)))+dx*dy/dt;
                u(i,j)=uu(i,j)-(pex(i+1,j)-pex(i,j))*dy/ac;
            end
        end
        for i=2:nx+2
             for j=2:ny+2
               bc=dx*(0.25*(v(i,j+1)-v(i,j-1))+2/(Re*dy))+dy*(0.25*(u(i,j+1)+u(i,j)-u(i-1,j+1)-u(i-1,j))+2/(Re*dx));
               v(i,j)=vv(i,j)-(pex(i,j+1)-pex(i,j))*dx/bc;
             end
        end
        
    %PPE    
    for i=2:nx+2
        for j=2:ny+2
           
            ac=dy*(0.25*(u(i+1,j)-u(i-1,j))+2/(Re*dx))+dx*(0.25*(v(i,j)+(v(i+1,j)-v(i,j-1)-v(i+1,j-1))+2/(Re*dy)))+dx*dy/dt;
            adown=-dy*(0.25*(u(i,j)+u(i-1,j))+1/(Re*dx));
            bc=dx*(0.25*(v(i,j+1)-v(i,j-1))+2/(Re*dy))+dy*(0.25*(u(i,j+1)+u(i,j)-u(i-1,j+1)-u(i-1,j))+2/(Re*dx));
            bleft=-dx*(0.25*(v(i,j-1)+v(i,j))+1/(Re*dy));
            
            cc=dy*dy/ac+dy*dy/aleft+dx*dx/bc+dx*dx/bleft;
            cup=dy*dy/ac;
            cdown=dy*dy/adown;
            cright=dx*dx/bc;
            cleft=dx*dx/bleft;
            ex(i,j)=-dy*(u(i,j)-u(i-1,j))-dx*(v(i,j)-v(i,j-1));
            
            ppex(i,j)=ex(i,j)+cleft*ppex(i,j-1)+cdown*ppex(i-1,j)+cup*pex(i+1,j)+cright*pex(i,j+1);
            
            
        end
    end
    
    pex=ppex;
    
    end
    val=max(max(abs(ex(2:end-1,2:end-1))));
    
    pp(2:end-1,2:end-1)=p(2:end-1,2:end-1)+pex(:,:);
    
    
end
%%
 quiver(x,y,u(2:end-1,2:end-1)',v(2:end-1,2:end-1)',0.6,'k-')
        axis equal
        axis([0 L 0 H])
        hold on 
        pcolor(x,y,p');
        colormap(jet)
        colorbar
        shading interp
        axis equal
        axis([0 L 0 H])
        title({['2D Cavity Flow with Re = ',num2str(Re)];['time(\itt) = ',num2str(dt*it)]})
        xlabel('Spatial co-ordinate (x) \rightarrow')
        ylabel('Spatial co-ordinate (y) \rightarrow')
        drawnow;
        hold off
