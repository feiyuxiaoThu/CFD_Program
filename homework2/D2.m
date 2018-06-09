%FTCS 格式求解二维问题
clear ;

d = 0.1; % 0.1,0.5,1.0

num_x = 200; % 100-500
num_y=num_x;
dx = 1.0/num_x;
dy = 1.0/num_y;
%dx=dy
dt = d * dx*dx;

t_total=0.01; %0.01,0.1,1,10

num_t=ceil(t_total/dt);
xx=zeros(1,num_x);
tt=zeros(1,num_t);
xx(1)=0;

for i = 2:1:num_x
xx(i) = xx(i-1) + dx;
end

tt(1)=0;
for i=2:1:num_t
tt(i) = tt(i-1) +dt;
end

dMx=0.01; %dMx=dx
dMt=0.001;
Mx=1.0/dMx; % 每过0.001记录一个值  100 
Mt=t_total/dMt;

M_u = zeros(num_x,num_y);%不一定要全部输入存储

dttnum=fix(dMt/dt); 


U1 = zeros(num_x,num_x);

U1(1,:)=0;
U1(:,1)=0;       
U1(num_x,:)=100;       
U1(:,num_y)=50;

U2 = zeros(num_x,num_x);

% 初始温度为0


  
        
for n = 2:1:num_t
    
        U1(1,:)=0;
        U1(:,1)=0;
        U1(num_x,:)=100;
        U1(:,num_y)=100;    %边界条件
        
        U2(:,:)=U1(:,:);
    
        for i=2:1:num_x-1
            for j=2:1:num_y-1
            
            U1(i,j)=d*U2(i+1,j) + (0.5-2*d)*U2(i,j) + d*U2(i-1,j) + d*U2(i,j+1) + (0.5-2*d)*U2(i,j) + d*U2(i,j-1) ;
            
            end
            
            
        end
        
        
        
       
        
      U2(:,:)=U1(:,:);
      
        
        
    
    
     
   
end

 M_u(:,:)=U2(:,:);

% plot(M_u);
 mesh(M_u);
title('T=0.01,Two Dimension');
xlabel('X');
ylabel('Y');



    

