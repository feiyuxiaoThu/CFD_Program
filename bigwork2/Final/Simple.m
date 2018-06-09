%Simulating a Lid-Driven Cavity flow with a lid moving at speed U
...Flow considered is laminar and of an incompressible Newtonian fluid
...Pressure correction method is used on a staggered grid
...to solve the Navier Stokes equations
%%
clear ;
clc;

%参数
L=1;         %方腔的边长
U=1;         %平板运动速度
Re=1000;      %雷诺数
dt=0.000001;     %时间步长
N=51;        %x方向的格点数
M=51;        %y方向的格点数
dx=1/(N-1);  %x方向的网格间距
dy=1/(M-1);  %y方向的网格间距


%交错网格下标参数的区分
i_u=1:(N-1);        %u in x direction
j_u=1:M;            %u in y direction
i_v=1:N;            %v in x direction
j_v=1:(M-1);        %v in y direction
i_p=1:N;            %p in x direction
j_p=1:M;            %p in y direction

i_uip=2:(N-2);    %u of inner points
j_uip=2:(M-1);      
i_vip=2:(N-1);    %v of inner points
j_vip=2:(M-2);
i_pip=2:(N-1);    %p of inner points
j_pip=2:(M-1);

%%
%内部格点赋初值
%由于为定常问题，初始值对于最终的结果没有影响，所以在满足方程的初始值中零值场为最简单的初始值
u(i_u,j_u)=0;         
v(i_v,j_v)=0;         
p(i_p,j_p)=0;         

% 给出初始估计值
u_star(1:(N-1),1:M)=0; u_star(:,:)=0;       
v_star(1:N,1:(M-1))=0; v_star(:,:)=0;       
p_star(1:N,1:M)=0; p_star(:,:)=0;           

%初始边界条件
u_star(i_u,M)=U;                 %无滑移条件
u_star(i_u,1)=0;                
v_star(1,j_v)=0;               
v_star(N,j_v)=0;                 

%Left Boundary velocity for u
u_star(1,j_u)=(u_star(2,j_u))/3;              %u velocity near wall approx to one third of adjacent u node

%Right Boundary velocity for u
u_star(N-1,j_u)=(u_star(N-2,j_u))/3;          %u velocity near wall approx to one third of adjacent u node
 
%Top Boundary velocity for v
v_star(i_v,M-1)=(v_star(i_v,M-2))/3;           %v velocity near top approx to one third of adjacent v node 
 
%Bottom Boundary velocity for v
v_star(i_v,1)=(v_star(i_v,2))/3;              %v velocity near top approx to one third of adjacent v node


%压力边界条件均采用向内零阶延拓的边界条件
p(i_p,1)=p(i_p,2);                  
p(1,j_p)=p(2,j_p);                  
p(N,j_p)=p(N-1,j_p);                
p(i_p,M)=p(i_p,(M-1));              

%%
%-------------------循环迭代直到收敛-----------------
% SIMPLE algorithm
%速度场和压力场的修正量      
u_corr(1:(N-1),1:M)=0; u_corr(:,:)=0;
v_corr(1:N,1:(M-1))=0; v_corr(:,:)=0;
p_corr(1:N,1:M)=0; p_corr(:,:)=0;

for i=2:N-1
    for j=2:M-1
        p_corr(i,j)=0.1;  %和后面的收敛条件对应即可
    end 
end


 while (max(max(abs(p_corr)))>0.00005)       % 收敛条件

% 在假设的压力分布下求解动量方程
	u_bar(1:(N-1),1:M)=0;
    for i=2:(N-1)
        for j=1:(M-1)
            u_bar(i,j)=0.25*(u_star(i,j)+u_star(i-1,j)+u_star(i,j+1)+u_star(i-1,j+1));      %Average values of u
        end
    end

    v_bar(1:N,1:(M-1))=0;
	for i=1:(N-1)
        for j=2:(M-1)
            v_bar(i,j)=0.25*(v_star(i,j)+v_star(i+1,j)+v_star(i,j-1)+v_star(i+1,j-1));      %Average values of v
        end
    end
    
    %x方向的动量方程
    u_new(1:(N-1),1:M)=0;
    A(1:(N-1),1:M)=0;
    for i=2:(N-2)
        for j=2:(M-1)
            A(i,j)=-(((u_star(i,j)*(u_star(i+1,j)-u_star(i-1,j)))/(2*dx))+((v_bar(i,j)*(u_star(i,j+1)-u_star(i,j-1)))/(2*dy)))+(1/Re)*(((u_star(i+1,j)-(2*u_star(i,j))+u_star(i-1,j))/(dx^2))+((u_star(i,j+1)-2*u_star(i,j)+u_star(i,j-1))/(dy^2)));
            u_new(i,j)=u_star(i,j)+(dt*A(i,j))-(dt*(p_star(i+1,j)-p_star(i,j)))/dx;
        end
    end
    
    u_new(1:(N-1),M)=U;                         
    u_new(1:(N-1),1)=0;                         
    for j=2:(M-1)                             
       u_new(1,j)=(u(2,j))/3;                 
      u_new(N-1,j)=(u(N-2,j))/3;

    end

     %y方向的动量方程
    v_new(1:N,1:(M-1))=0;
    B(1:N,1:(M-1))=0;
    for i=2:(N-1)
        for j=2:(M-2)
            B(i,j)=-(((v_star(i,j)*(v_star(i,j+1)-v_star(i,j-1)))/(2*dy))+((u_bar(i,j)*(v_star(i+1,j)-v_star(i-1,j)))/(2*dx)))+(1/Re)*(((v_star(i+1,j)-(2*v_star(i,j))+v_star(i-1,j))/(dx^2))+((v_star(i,j+1)-2*v_star(i,j)+v_star(i,j-1))/(dy^2)));
            v_new(i,j)=v_star(i,j)+(dt*B(i,j))-(dt*(p_star(i,j+1)-p_star(i,j)))/dy;
        end
    end

    v_new(1,1:(M-1))=0;                         
    v_new(N,1:(M-1))=0;                         
    for i=2:(N-1)                               
         v_new(i,M-1)=(v_new(i,M-2))/3;          
        v_new(i,1)=(v_new(i,2))/3;     

    end
    
% 求解压力修正方程
% 迭代法求解泊松方程
        
    p_cor1(1:N,1:M)=0;
    for r=1:20   %迭代次数
        for i=2:(N-1)
            for j=2:(M-1)
                D(i,j)=(1/dt)*(((u_new(i,j)-u_new(i-1,j))/dx)+((v_new(i,j)-v_new(i,j-1))/dy));
                p_cor1(i,j)=(((dx^2)*(dy^2))/(2*(dy^2)+2*(dx^2)))*(((p_corr(i+1,j)+p_corr(i-1,j))/(dx^2))+((p_corr(i,j+1)+p_corr(i,j-1))/(dy^2))-D(i,j));
            end
        end
        p_corr=p_cor1;
    end

% 得到修正后的压力场        

    for i=1:N
        for j=1:M
            p(i,j)=p_star(i,j)+p_corr(i,j);
        end
    end

    %边界条件的修正
    
    p(i_p,1)=p(i_p,2);                  
    p(1,j_p)=p(2,j_p);                  
    p(N,j_p)=p(N-1,j_p);                
    p(i_p,M)=p(i_p,(M-1));              
    
    p_star=p;
    
% 利用修正的压力场求解修正速度场     
    
    for i=1:(N-1)
        for j=1:M
            u_corr(i,j)=-dt*((p_corr(i+1,j)-p_corr(i,j))/dx);
            u(i,j)=u_star(i,j)+u_corr(i,j);
            
            
            u(1:(N-1),M)=U;                        
            u(1:(N-1),1)=0;                         

            for j_u=2:(M-1)                       
                 u(1,j_u)=(u(2,j_u))/3;                                
            end

            for j_u=2:(M-1)                       
                u(N-1,j_u)=(u(N-2,j_u))/3;                      
            end
        end
    end
    
    u_star=u;
    
    for i=1:N
        for j=1:(M-1)
            v_corr(i,j)=-dt*((p_corr(i,j+1)-p_corr(i,j))/dy);
            v(i,j)=v_star(i,j)+v_corr(i,j);
            
            
             v(1,1:(M-1))=0;                         
             v(N,1:(M-1))=0;                       

            for i_v=2:(N-1)                         
                v(i_v,M-1)=(v(i_v,M-2))/3;       
            end

            for i_v=2:(N-1)                        
                v(i_v,1)=(v(i_v,2))/3;               
            end
        end
    end
       
    v_star=v;

end

%------------------循环结束--------------------

%%
% 计算非结构网格上的速度压力场
%u_plot(1:N,1:M)=0; u_plot(:,:)=0;

for i=2:N-1
    for j=2:M-1
        u_plot(i,j)=0.5*(u_star(i-1,j)+u_star(i,j));
    end
end

u_plot(2:(N-1),M)=U;                
u_plot(1:(N-1),1)=0;                
u_plot(1,1:M)=0;                   
u_plot(N,1:M)=0;                    
u_plot=u_plot';

for i=2:N-1
    for j=2:M-1
        v_plot(i,j)=0.5*(v_star(i,j-1)+v_star(i,j));
    end
end

v_plot(2:(N-1),M)=0;                
v_plot(1:(N-1),1)=0;                
v_plot(1,1:M)=0;                    
v_plot(N,1:M)=0;                    
v_plot=v_plot';

for i=1:N
    for j=1:M
        p_plot(i,j)=p(i,j);
    end
end
%%
x=zeros(M,N);
xx=zeros(M);
y=zeros(M,N);
yy=zeros(N);
for i=1:M
    x(i,:)=dx*(i-1);
    xx(i)=dx*(i-1);
end
for j=1:N
    y(:,j)=dy*(j-1);
    yy(j)=dy*(j-1);
end


quiver(u_plot,v_plot,0.6,'k-');

        pcolor(p);
        colormap(jet)
        colorbar
        shading interp
        xlabel('Spatial co-ordinate (x) \rightarrow')
        ylabel('Spatial co-ordinate (y) \rightarrow')
        drawnow;
        hold off
        
        
figure
streamslice(u_plot,v_plot);

figure
mesh(p);