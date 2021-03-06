%FTCS 格式
clear ;

d = 0.1; % 0.1,0.5,1.0

num_x = 200; % 100-500
dx = 1.0/num_x;
dt = d * dx*dx;

t_total=10; %0.01,0.1,1,10

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
dMt=0.0001;
Mx=1.0/dMx; % 每过0.001记录一个值  100 
Mt=t_total/dMt;
M_u = zeros(num_x,Mt);%不一定要全部输入存储

dttnum=fix(dMt/dt); 


U1x = zeros(1,num_x);
U2x = zeros(1,num_x);


 for i = 2:1:num_x-1
            %U1x(i) = fx(xx(i));
            
            if xx(i)<0.3
                ut0 = 0;
            elseif  (xx(i)<=0.6) && (xx(i)>=0.3)
             ut0 = 1;
            else
             ut0 = 1 + 2.5* (xx(i) - 0.6);
            end
            
           
           U1x(i)=ut0;
            U2x(i) = U1x(i);
            
 end
        
         U1x(1)=0;
                 
        U1x(num_x)=bt(2,tt(1)); % 1,2 choice
        U2x(1)=U1x(1);
        U2x(num_x)=U1x(num_x);
        
        M_u(:,1)=U2x(:);
        Mii=1;
        
for n = 2:1:num_t
    
    
    
        for i=2:1:num_x-1
            
            U1x(i)=d*U2x(i+1) + (1-2*d)*U2x(i) + d*U2x(i-1);
            
        end
        
        
        
        U1x(1)=0;
        
        U1x(num_x)=bt(2,tt(n)); % 1,2 choice
        
      
        U2x(:)=U1x(:);
        
      
        
        
     if (mod(n,dttnum)==0)
        Mii=Mii+1;
      M_u(:,Mii)=U2x(:);
     end
    
   
end


% plot(M_u);
 mesh(M_u);
title('T=10,With boundary condition2 and FTCS');
xlabel('T');
ylabel('X');



    

