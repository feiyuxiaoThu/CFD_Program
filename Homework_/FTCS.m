%FTCS 格式
clear ;

d = 0.5; % 0.1,0.5,1.0

num_x = 500; % 100-500
dx = 1.0/num_x;
dt = d * dx*dx;
num_t = 100;
t_total=ceil(num_t*dt);

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


M_u = zeros(num_x,num_t);%不一定要全部输入存储

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
                 
        U1x(num_x)=bt(1,tt(1)); % 1,2 choice
        U2x(1)=U1x(1);
        U2x(num_x)=U1x(num_x);
        
        M_u(:,1)=U2x(:);
        
        
for n = 2:1:num_t
    
    
    
        for i=2:1:num_x-1
            
            U1x(i)=d*U2x(i+1) + (1-2*d)*U2x(i) + d*U2x(i-1);
            
        end
        
        
        
        U1x(1)=0;
        
        U1x(num_x)=bt(1,tt(n)); % 1,2 choice
        
      
        U2x(:)=U1x(:);
        
        M_u(:,n)=U2x(:);
        
        
    
    
   
end


mesh(M_u);
title('The Temperature distribution as time evolves,dt = ');
xlabel('T');
ylabel('X');



    

