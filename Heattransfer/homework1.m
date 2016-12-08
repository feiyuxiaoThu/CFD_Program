%传热学作业 
%FTCS 格式
clear ;


a=9/17600;

dt = 3600;

t_total=24;
dx=0.01;

d=dt/(dx*dx*a);

x_total=0.24;
num_x=x_total/dx+1;

num_t=t_total;
xx=zeros(1,num_x);
tt=zeros(1,num_t);

%初始化 
xx(1)=0;

for i = 2:1:num_x
xx(i) = xx(i-1) + dx;
end

tt(1)=0;
for i=2:1:num_t
tt(i) = tt(i-1) +dt;
end




M1_u = zeros(num_x,3);
M2_u = zeros(num_x,3);%对比前后




U1x = zeros(1,num_x);
for i=1:num_x
    U1x(i)=15;
end
U2x = zeros(1,num_x);
U2x(:)=U1x(:);

U0=[-5.9,-6.2,-6.6,-6.7,-6.8,-6.9,-7.2,-7.7,-7.6,-7.0,-4.9,-2.3,-1.0,-2.4,1.8,1.8,1.6,0.5,-1.6,-2.8,-3.5,-4.3,4.8,-5.3];
 
    
    
      
        
  h1=10;
  h2=6;
  nad=0.81;
  b1=h1*dx/nad;
  b2=h2*dx/nad;
  
  
  M=0;
    exit=0;
    
  while(exit==0)
      
      
    
  
    for n=1:1:24  %控制时间
        
     %周期边界条件
     Uf2=15;
     ii=n;
     Uf1=U0(ii);
     
     U1x(1)=(2*d*(1+b1)+1)*U2x(1)-2*d*U2x(2)-2*d*b1*Uf1;
     U1x(num_x)=(2*d*(1+b1)+1)*U2x(num_x)-2*d*U2x(num_x-1)-2*d*b1*Uf2;
     
      M1_u(ii,1)=U1x(1);
      M1_u(ii,2)=U1x((1+num_x)/2);
      M1_u(ii,3)=U1x(num_x);
        
      M2_u(:,:)=M1_u(:,:);
                           
                        
        for i=2:1:num_x-1
            
            U1x(i)=d*U2x(i+1) + (1-2*d)*U2x(i) + d*U2x(i-1);
            
        end
       
%         M1_u = zeros(3,num_x);
%         M2_u = zeros(3,num_x);%对比前后

        M2_u(ii,1)=U1x(1);
        M2_u(ii,2)=U1x((1+num_x)/2);
        M2_u(ii,3)=U1x(num_x);
        
      
          U2x(:)=U1x(:); 
        
        
        
        
    end
    
    
    
    
 for i1=1:3
        for i2=1:num_x
          M=abs(M1_u(i2,i1)-M2_u(i2,i1)) + M;
        end
  end
    
    if M<0.1
       exit=1;
    end
       
    
    
  end
  
     
    
   







    



