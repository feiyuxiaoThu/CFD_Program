%% 说明
%Simulating a Lid-Driven Cavity flow with a lid moving at speed U_move
...Flow considered is laminar and of an incompressible Newtonian fluid
...Pressure correction method is used on a staggered grid
...to solve the Navier Stokes equations

clear
clc
close all

%% 计算参数和变量
rho=1; %流体密度（视为常数）
mu=1e-5; % 水的粘性系数（取其数量级）
% 计算区域
M=20; % x方向网格数
N=M;  %正方形区域
L=0.001;
H=L;
dx=L/M; % 网格宽度
dy=dx; % 方形网格
% Convergence ？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
alpha_p=.001; % Pressure under-relaxtion
alpha_u=.02; % u-momentum under-relaxation
alpha_v=.02; % v-momentum under-relaxation
max_iter=20000; % 最大迭代次数
max_res=1e-4; % 收敛判据
inner_iter=1; % Number of iterations before velocity is updated
% Parameters
Dn=mu*dx/dy;Ds=Dn;
De=mu*dy/dx;Dw=De;
%雷诺数
Re=4000; % 400,1000

% Boundary Conditions
P_free=10; % 压强初始给定一个均匀场实际最后求得的压强是一个相对值与初始给定的压强无关
u_move=Re*mu/L;


%% 初始速度场赋值
% 初始值只需要满足边界条件和连续方程即可，对最终收敛的结果也不会产生影响
u_ini=zeros(M+2,N+1); % m/s
u_ini(1,:)=u_move*ones(1,N+1);
u_ini(2:end-1,2:end-1)=u_move*ones(size(u_ini)-2);
u_ini=(rot90(u_ini,3));

v_ini=zeros(M+1,N+2); 
v_ini(2:end-1,2:end-1)=u_move*ones(size(v_ini)-2);
v_ini=(rot90(v_ini,3)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 迭代前值
u_old=u_ini; 
v_old=v_ini; 
% 迭代后值
u_new=u_ini; % m/s
apm_u=zeros(size(u_new)); % The u momentum :a" coefficient matrix
v_new=v_ini; % m/s
apm_v=zeros(size(v_new)); % The v momentum "a" coefficient matrix 
% 求压强
x=linspace(0,1,M);
y=x;
[x , y]=meshgrid(x,y);
z=1./((x-1).^2+(y).^2+.5);
P=P_free*z; % Pascals
P=(rot90(P,3)); % Pascals, rotating to get i corresponding 

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
    count=count+1; 
    disp(count)
    for k=1:inner_iter
    %% 在假定压力场下估算速度
    for i=2:size(u_new,1)-1
        for j=2:size(u_new,2)-1
            
            % Computing the face fluxes:
            Fe=dy/2*(u_old(i+1,j)+u_old(i,j)); Fw=dy/2*(u_old(i,j)+u_old(i-1,j));
            Fn=dx/2*(v_old(i,j)+v_old(i+1,j)); Fs=dx/2*(v_old(i,j-1)+v_old(i+1,j-1));
            
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
    for i=2:size(v_new,1)-1
        for j=2:size(v_new,2)-1
            
            % Computing the face fluxes:
            Fe=dy/2*(u_old(i,j+1)+u_old(i,j)); Fw=dy/2*(u_old(i-1,j+1)+u_old(i-1,j));
            Fn=dx/2*(v_old(i,j+1)+v_old(i,j)); Fs=dx/2*(v_old(i,j)+v_old(i,j-1));
            
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
%     v_s(2:end-1,end)=v_s(2:end-1,end-1);
    %     for i=2:N+1
%         u_s(i,M+2)=u_s(i-1,M+2)+v_s(i,M+1);
%     end


    %% 求解压力修正方程

    % Pressure correction loop
    for I=1:N
        for J=1:M
           if I==1 && J==1
               % Bottom left corner
               ae=dy^2/apm_u(I+1,J+1); aw=0;
               an=dx^2/apm_v(I+1,J+1); as=0;
               ap=ae+aw+an+as;
               b=dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1);
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+an*P_p(I,J+1)+b);
           elseif I==1 && J~=1 && J~=M
               % Left wall
               ae=dy^2/apm_u(I+1,J+1); aw=0;
               an=dx^2/apm_v(I+1,J+1); as=dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1);
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+an*P_p(I,J+1)+as*P_p(I,J-1)+b);
           elseif J==1 && I~=1 && I~=N
               % Bottom wall
               ae=dy^2/apm_u(I+1,J+1); aw=dy^2/apm_u(I,J+1);
               an=dx^2/apm_v(I+1,J+1); as=0;
               ap=ae+aw+an+as;
               b=dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1);
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+aw*P_p(I-1,J)+an*P_p(I,J+1)+b);
           elseif I==N && J==1
               % Bottom right corner
               ae=0; aw=dy^2/apm_u(I,J+1);
               an=dx^2/apm_v(I+1,J+1); as=0;
               ap=ae+aw+an+as;
               b=dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1);
               P_p(I,J)=1/ap*(aw*P_p(I-1,J)+an*P_p(I,J+1)+b);
           elseif I==N && J~=1 && J~=M
               % Right wall
               ae=0; aw=dy^2/apm_u(I,J+1);
               an=dx^2/apm_v(I+1,J+1); as=dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1);
               P_p(I,J)=1/ap*(aw*P_p(I-1,J)+an*P_p(I,J+1)+as*P_p(I,J-1)+b);
           elseif I==1 && J==M
               % Top Left
               ae=dy^2/apm_u(I+1,J+1); aw=0;
               an=0; as=dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1);
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+as*P_p(I,J-1)+b);
           elseif I==N && J==M
               % Top Right
               ae=0; aw=dy^2/apm_u(I,J+1);
               an=0; as=dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1);
               P_p(I,J)=1/ap*(aw*P_p(I-1,J)+as*P_p(I,J-1)+b);
           elseif I~=N && I~=1 && J==M
               % Top Row
               ae=dy^2/apm_u(I+1,J+1); aw=dy^2/apm_u(I,J+1);
               an=0; as=dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=dy*u_ini(I,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1);
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+aw*P_p(I-1,J)+as*P_p(I,J-1)+b);
           else
               % Interior nodes
               ae=dy^2/apm_u(I+1,J+1); aw=dy^2/apm_u(I,J+1);
               an=dx^2/apm_v(I+1,J+1); as=dx^2/apm_v(I+1,J);
               ap=ae+aw+an+as;
               b=dy*u_ini(I,J+1)-dy*u_ini(I+1,J+1)+dx*v_ini(I+1,J)-dx*v_ini(I+1,J+1);
               P_p(I,J)=1/ap*(ae*P_p(I+1,J)+aw*P_p(I-1,J)+an*P_p(I,J+1)+as*P_p(I,J-1)+b);
                
            end % the conditional boundary node statement
        end % the J loop
    end % the I loop

    %% 修正压力场合速度场
    P=P_s+alpha_p*P_p; % Relaxed Pressure
    
    % u-correction
    for i=2:size(u_new,1)-1
        for j=2:size(u_new,2)-1
            u_new(i,j)=u_ini(i,j)+dy/apm_u(i,j)*(P_p(i-1,j-1)-P_p(i,j-1));
        end
    end
    % v-correction
    for i=2:size(v_new,1)-1
        for j=2:size(v_new,2)-1
            v_new(i,j)=v_ini(i,j)+dx/apm_v(i,j)*(P_p(i-1,j-1)-P_p(i-1,j));
        end
    end
    %% Enforcing conservation of mass on the top boundary
%     v(2:end-1,end)=v(2:end-1,end-1);

    %% Relaxing the velocities
    % Similar equations to "solving the starred momentum", only the starred
    % quantities are the new guesses (both velocity and pressure) 
    
    % u-momentum
    for i=2:size(u_new,1)-1
        for j=2:size(u_new,2)-1
            
            % Computing the face fluxes:
            Fe=dy/2*(u_old(i+1,j)+u_old(i,j)); Fw=dy/2*(u_old(i,j)+u_old(i-1,j));
            Fn=dx/2*(v_old(i,j)+v_old(i+1,j)); Fs=dx/2*(v_old(i,j-1)+v_old(i+1,j-1));
            
            % Computing the "a" coefficients:
            ae=De-min(0,Fe); aw=Dw+max(0,Fw);
            an=Dn-min(0,Fn); as=Ds+max(0,Fs);
            ap=ae+aw+an+as+(Fe-Fw+Fn-Fs);
            
            u_new(i,j)=alpha_u/ap*(ae*u_new(i+1,j)+aw*u_new(i-1,j)... % x-direction
                +an*u_new(i,j+1)+as*u_new(i,j-1)...% y-direction
                -dy*(P(i,j-1)-P(i-1,j-1)))+(1-alpha_u)*u_old(i,j); % Pressure term
            
        end % the j loop
    end % the i loop
    
    % v-momentum
    for i=2:size(v_new,1)-1
        for j=2:size(v_new,2)-1

            % Computing the face fluxes:
            Fe=dy/2*(u_old(i,j+1)+u_old(i,j)); Fw=dy/2*(u_old(i-1,j+1)+u_old(i-1,j));
            Fn=dx/2*(v_old(i,j+1)+v_old(i,j)); Fs=dx/2*(v_old(i,j)+v_old(i,j-1));

            % Computing the "a" coefficients:
            ae=De-min(0,Fe); aw=Dw+max(0,Fw);
            an=Dn-min(0,Fn); as=Ds+max(0,Fs);
            ap=ae+aw+an+as+(Fe-Fw+Fn-Fs);

            v_new(i,j)=alpha_v/ap*(ae*v_new(i+1,j)+aw*v_new(i-1,j)... % x-direction
                +an*v_new(i,j+1)+as*v_new(i,j-1)...% y-direction
                -dx*(P(i-1,j)-P(i-1,j-1)))+(1-alpha_v)*v_old(i,j); % Pressure term


        end % the j loop
    end % the I loop
     %% Enforcing conservation of mass on the top boundary
%         v(2:end-1,end)=v(2:end-1,end-1);
    %     for i=2:N+1
%         u(i,M+2)=u(i-1,M+2)+v(i,M+1);
%     end
    end % the inner iterations
     %% 结束迭代
    if count>max_iter
        disp('Maximum number of iterations has been reached')
        break 
    end
    if max(abs(u_old(:)-u_new(:)))<max_res && max(abs(v_old(:)-v_new(:)))<max_res && max(abs((P_o(:)-P(:))./P_o(:)))<10*max_res
        disp('Success: Solution has converged');
        break  
    end
    %% 准备下次迭代赋值
    u_ini=u_new;
    v_ini=v_new;
    P_s=P;
    u_old=u_new;
    v_old=v_new;
    P_o=P;
end

%% 数据可视化（绘制速度和压力分布图）
x=linspace(0,1,N+1);
y=x;
[x , y]=meshgrid(x,y);
figure;quiver(x,y,u_new(:,2:end)',v_new(2:end,:)'',0.6,'k-'); hold on ; colormap(jet);colorbar;shading interp;xlabel('Spatial co-ordinate (x) \rightarrow')
% figure,quiver(x,y,u(:,2:end)',v(2:end,:)',2,'k');hold on
% xlabel('x','fontsize',16);
% ylabel('y','fontsize',16);
% set(gca,'fontsize',16);

x=linspace(0,1,N);
y=x;
[x , y]=meshgrid(x,y);
figure;
pcolor(x,y,P');shading interp
pco=get(gca,'children');
set(pco(1),'facealpha',.5);
xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
set(gca,'fontsize',16);
colorbar
% figure,surf(u);shading interp
% figure,surf(v);shading interp

 