%Simulating a Lid-Driven Cavity flow with a lid moving at speed U_move
...Flow considered is laminar and of an incompressible Newtonian fluid
...Pressure correction method is used on a staggered grid
...to solve the Navier Stokes equations

clear
clc
close all

%% 计算参数和常量
rho=1; % 流体（水）的密度
mu=1e-5; % 粘性系数（取均值——数量级）
% 计算区域
nx=50; % x方向网格数
ny=nx;
L=1;
H=L;
dx=L/nx; 
dy=H/ny; 
% Convergence
alpha_p=.001; % Pressure under-relaxtion
alpha_u=.02; % u-momentum under-relaxation
alpha_v=.02; % v-momentum under-relaxation
max_iter=200000; % 最大迭代次数
max_res=1e-6; % 收敛判据
inner_iter=1; % Number of iterations before velocity is updated
% Parameters
Dn=mu*dx/dy;Ds=Dn;
De=mu*dy/dx;Dw=De;

%雷诺数
Re=400; % 400,1000
U_move=Re*mu/L;


%% 边界条件
P_free=10; % 压强初始给定一个均匀场实际最后求得的压强是一个相对值与初始给定的压强无关
u_free=U_move; % meters/second


%% 初始速度场赋值
% 初始值只需要满足边界条件和连续方程即可，对最终收敛的结果也不会产生影响
u_ini=zeros(nx+2,ny+1); 
u_ini(1,:)=u_free*ones(1,ny+1);
u_ini(2:end-1,2:end-1)=u_free*ones(size(u_ini)-2);
u_ini=(rot90(u_ini,3)); 
v_ini=zeros(nx+1,ny+2); 
v_ini(2:end-1,2:end-1)=u_free*ones(size(v_ini)-2);
v_ini=(rot90(v_ini,3)); 



% 迭代前
u=u_ini; 
v=v_ini; 
% 迭代后
uu=u_ini; 
apm_u=zeros(size(uu)); % The u momentum :a" coefficient matrix
vv=v_ini; 
apm_v=zeros(size(vv)); % The v momentum "a" coefficient matrix 
% Acutal Pressure
x=linspace(0,1,nx);
y=linspace(0,1,ny);
[x , y]=meshgrid(x,y);
z=1./((x-1).^2+(y).^2+.5);
P=P_free*z; % Pascals
P=(rot90(P,3)); 

% Starred Pressure

P_s=P; % Pascals
% Old Pressure
P_o=P; % Pascals
% Pressure Correction
P_p=zeros(size(P)); % Pascals

%% Initiating the Simple Algorithm Infinite Loop
% Begin the while loop
count=0;
residual=zeros(3,max_iter);

while 1
    count=count+1; % 计数器
    disp(count)
    for k=1:inner_iter
    %% Solving starred momentum equations
%     u=u+.5*exp(-0.05*count);
%     v=v+5*exp(-0.1*count);    
%     P=P+50*exp(-count);
    % u-momentum
    for i=2:size(uu,1)-1
        for j=2:size(uu,2)-1
            
            % Computing the face fluxes:
            Fe=rho*dy/2*(u(i+1,j)+u(i,j)); Fw=rho*dy/2*(u(i,j)+u(i-1,j));
            Fn=rho*dx/2*(v(i,j)+v(i+1,j)); Fs=rho*dx/2*(v(i,j-1)+v(i+1,j-1));
            
            % Computing the "a" coefficients:
            ae=De-min(0,Fe); aw=Dw+max(0,Fw);
            an=Dn-min(0,Fn); as=Ds+max(0,Fs);
            apm_u(i,j)=ae+aw+an+as+(Fe-Fw+Fn-Fs);
            
            u_ini(i,j)=1/apm_u(i,j)*(ae*u_ini(i+1,j)+aw*u_ini(i-1,j)... % x-direction
                +an*u_ini(i,j+1)+as*u_ini(i,j-1)...% y-direction
                -dy*(P_s(i,j-1)-P_s(i-1,j-1))); % Pressure term
            
        end % the j loop
    end % the i loop
    
    % v-momentum
    for i=2:size(vv,1)-1
        for j=2:size(vv,2)-1
            
            % Computing the face fluxes:
            Fe=rho*dy/2*(u(i,j+1)+u(i,j)); Fw=rho*dy/2*(u(i-1,j+1)+u(i-1,j));
            Fn=rho*dx/2*(v(i,j+1)+v(i,j)); Fs=rho*dx/2*(v(i,j)+v(i,j-1));
            
            % Computing the "a" coefficients:
            ae=De-min(0,Fe); aw=Dw+max(0,Fw);
            an=Dn-min(0,Fn); as=Ds+max(0,Fs);
            apm_v(i,j)=ae+aw+an+as+(Fe-Fw+Fn-Fs);
            
            v_ini(i,j)=1/apm_v(i,j)*(ae*v_ini(i+1,j)+aw*v_ini(i-1,j)... % x-direction
                +an*v_ini(i,j+1)+as*v_ini(i,j-1)...% y-direction
                -dx*(P_s(i-1,j)-P_s(i-1,j-1))); % Pressure term
            
            
        end % the j loop
    end % the I loop
    
    %% Enforcing conservation of mass on the top boundary
     v_ini(2:end-1,end)=v_ini(2:end-1,end-1);
         for i=2:nx+1
                u_ini(i,nx+2)=u_ini(i-1,nx+2)+v_ini(i,nx+1);
        end


    %% Solving the pressure correction equations

    % Pressure correction loop
    for I=1:ny
        for J=1:nx
           if I==1 && J==1
               % Bottom left corner
               ae=rho*dy^2/apm_u(I+1,J+1); aw=0;
               an=rho*dx^2/apm_v(I+1,J+1); as=0;
               ap=ae+aw+an+as;
               b=rho*(dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1));
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+an*P_p(I,J+1)+b);
           elseif I==1 && J~=1 && J~=nx
               % Left wall
               ae=rho*dy^2/apm_u(I+1,J+1); aw=0;
               an=rho*dx^2/apm_v(I+1,J+1); as=rho*dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=rho*(dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1));
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+an*P_p(I,J+1)+as*P_p(I,J-1)+b);
           elseif J==1 && I~=1 && I~=ny
               % Bottom wall
               ae=rho*dy^2/apm_u(I+1,J+1); aw=rho*dy^2/apm_u(I,J+1);
               an=rho*dx^2/apm_v(I+1,J+1); as=0;
               ap=ae+aw+an+as;
               b=rho*(dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1));
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+aw*P_p(I-1,J)+an*P_p(I,J+1)+b);
           elseif I==ny && J==1
               % Bottom right corner
               ae=0; aw=rho*dy^2/apm_u(I,J+1);
               an=rho*dx^2/apm_v(I+1,J+1); as=0;
               ap=ae+aw+an+as;
               b=rho*(dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1));
               P_p(I,J)=1/ap*(aw*P_p(I-1,J)+an*P_p(I,J+1)+b);
           elseif I==ny && J~=1 && J~=nx
               % Right wall
               ae=0; aw=rho*dy^2/apm_u(I,J+1);
               an=rho*dx^2/apm_v(I+1,J+1); as=rho*dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=rho*(dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1));
               P_p(I,J)=1/ap*(aw*P_p(I-1,J)+an*P_p(I,J+1)+as*P_p(I,J-1)+b);
           elseif I==1 && J==nx
               % Top Left
               ae=rho*dy^2/apm_u(I+1,J+1); aw=0;
               an=0; as=rho*dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=rho*(dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1));
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+as*P_p(I,J-1)+b);
           elseif I==ny && J==nx
               % Top Right
               ae=0; aw=rho*dy^2/apm_u(I,J+1);
               an=0; as=rho*dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=rho*(dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1));
               P_p(I,J)=1/ap*(aw*P_p(I-1,J)+as*P_p(I,J-1)+b);
           elseif I~=ny && I~=1 && J==nx
               % Top Row
               ae=rho*dy^2/apm_u(I+1,J+1); aw=rho*dy^2/apm_u(I,J+1);
               an=0; as=rho*dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=rho*(dy*u_ini(I,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1));
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+aw*P_p(I-1,J)+as*P_p(I,J-1)+b);
           else
               % Interior nodes
               ae=rho*dy^2/apm_u(I+1,J+1); aw=rho*dy^2/apm_u(I,J+1);
               an=rho*dx^2/apm_v(I+1,J+1); as=rho*dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=rho*(dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1));
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+aw*P_p(I-1,J)+an*P_p(I,J+1)+as*P_p(I,J-1)+b);
                
            end % the conditional boundary node statement
        end % the J loop
    end % the I loop

    %% Correcting the pressure and velocities
    P=P_s+alpha_p*P_p; % Relaxed Pressure
    
    % u-correction
    for i=2:size(uu,1)-1
        for j=2:size(uu,2)-1
            uu(i,j)=u_ini(i,j)+dy/apm_u(i,j)*(P_p(i-1,j-1)-P_p(i,j-1));
        end
    end
    % v-correction
    for i=2:size(vv,1)-1
        for j=2:size(vv,2)-1
            vv(i,j)=v_ini(i,j)+dx/apm_v(i,j)*(P_p(i-1,j-1)-P_p(i-1,j));
        end
    end
    %% Enforcing conservation of mass on the top boundary
    v(2:end-1,end)=v(2:end-1,end-1);

    %% Relaxing the velocities
    % Similar equations to "solving the starred momentum", only the starred
    % quantities are the new guesses (both velocity and pressure) 
    
    % u-momentum
    for i=2:size(uu,1)-1
        for j=2:size(uu,2)-1
            
            % Computing the face fluxes:
            Fe=rho*dy/2*(u(i+1,j)+u(i,j)); Fw=rho*dy/2*(u(i,j)+u(i-1,j));
            Fn=rho*dx/2*(v(i,j)+v(i+1,j)); Fs=rho*dx/2*(v(i,j-1)+v(i+1,j-1));
            
            % Computing the "a" coefficients:
            ae=De-min(0,Fe); aw=Dw+max(0,Fw);
            an=Dn-min(0,Fn); as=Ds+max(0,Fs);
            ap=ae+aw+an+as+(Fe-Fw+Fn-Fs);
            
            uu(i,j)=alpha_u/ap*(ae*uu(i+1,j)+aw*uu(i-1,j)... % x-direction
                +an*uu(i,j+1)+as*uu(i,j-1)...% y-direction
                -dy*(P(i,j-1)-P(i-1,j-1)))+(1-alpha_u)*u(i,j); % Pressure term
            
        end % the j loop
    end % the i loop
    
    % v-momentum
    for i=2:size(vv,1)-1
        for j=2:size(vv,2)-1

            % Computing the face fluxes:
            Fe=rho*dy/2*(u(i,j+1)+u(i,j)); Fw=rho*dy/2*(u(i-1,j+1)+u(i-1,j));
            Fn=rho*dx/2*(v(i,j+1)+v(i,j)); Fs=rho*dx/2*(v(i,j)+v(i,j-1));

            % Computing the "a" coefficients:
            ae=De-min(0,Fe); aw=Dw+max(0,Fw);
            an=Dn-min(0,Fn); as=Ds+max(0,Fs);
            ap=ae+aw+an+as+(Fe-Fw+Fn-Fs);

            vv(i,j)=alpha_v/ap*(ae*vv(i+1,j)+aw*vv(i-1,j)... % x-direction
                +an*vv(i,j+1)+as*vv(i,j-1)...% y-direction
                -dx*(P(i-1,j)-P(i-1,j-1)))+(1-alpha_v)*v(i,j); % Pressure term


        end % the j loop
    end % the I loop
     %% Enforcing conservation of mass on the top boundary
%         v(2:end-1,end)=v(2:end-1,end-1);
    %     for i=2:N+1
%         u(i,M+2)=u(i-1,M+2)+v(i,M+1);
%     end
    end % the inner iterations
    %% Checking for convergence
    if sum(isnan(uu(:)))+sum(isnan(vv(:)))+sum(isnan(P(:)))>=1
        disp('Error: Solution has diverged');
        break % End the simple algorithm
    end
    if count>max_iter
        disp('Maximum number of iterations has been reached')
        break % End the simple algorithm
    end
    if max(abs(u(:)-uu(:)))<max_res && max(abs(v(:)-vv(:)))<max_res && max(abs((P_o(:)-P(:))./P_o(:)))<10*max_res
        disp('Success: Solution has converged');
        break % End the simple algorithm
    end
    %% Re-initiating old values
    % Since the algorithm did not converge, the new values are now the old
    % values
    u_ini=uu;
    v_ini=vv;
    P_s=P;
    u=uu;
    v=vv;
    P_o=P;
end

%% Plotting Results
% x=linspace(0,1,ny+1);
% y=x;
% [x , y]=meshgrid(x,y);
% 
% quiver(uu(2:end,2:end-1)',vv(2:end-1,2:end)',0.6,'k-')
%         axis equal
%         
%         hold on 
%         pcolor(x,y,P');
%         colormap(jet)
%         colorbar
%         shading interp
%        
%         
%     
%         xlabel('Spatial co-ordinate (x) \rightarrow')
%         ylabel('Spatial co-ordinate (y) \rightarrow')
%         drawnow;
%         hold off
