%二维槽道流动的边值问题
D=1;
dy=0.02; %标准
L=5;
dx=dy;
mu=10^-6; %粘性系数

%dt=dx*dx/5;
dt=0.001;
t=1;
t_total=t/dt;

Re=200;
U=Re*mu/L ;

num_y=D/dy+1;
num_x=L/dx+1;


y=0:dy:D;


maxa=0;

W=zeros(num_x,num_y);
S=zeros(num_x,num_y);
u=zeros(num_x,num_y);
v=zeros(num_x,num_y);
dss=zeros(num_x,num_y);

%边界条件
%壁面
S(:,1)=0;
S(:,num_y)=U*D;
u(:,1)=0;
u(:,num_y)=0;
v(:,1)=0;
v(:,num_y)=0;

for i=1:num_x
W(i,1)=2*S(i,1)/(dy*dy);
W(i,num_y)=(2*S(1,num_y-1)-2*U*D)/(dy*dy);
end

%进口
for j=2:num_y-1
    S(1,j)=U*dy*(i-1);
    u(1,j)=U;
end

%u(1,1)=U;
%u(1,num_y)=U;

for j=2:num_y-1
    W(1,j)=2*(S(2,j)-S(1,j))/(dx*dx)+(S(1,j+1)-2*S(1,j)+S(1,j-1))/(dy*dy);
end



for t=1:dt:t_total
    


WW=W;
SS=S;



%出口
for j=2:num_y-1
    SS(num_x,j)=2*S(num_x-1,j)-S(num_x-2,j);
end

for j=2:num_y-1
    WW(num_x,j)=W(num_x-1,j);
end



e=0.1;

%FTCS格式求解涡量方程
for i=2:num_x-1
    for j=2:num_y-2
        W(i,j)=( WW(i+1,j)-2*WW(i,j)+WW(i-1,j)+WW(i,j+1)-2*WW(i,j)+WW(i,j-1) )*dt/(Re*dx*dx)+WW(i,j)-( u(i,j)*( WW(i+1,j)-WW(i-1,j) )+v(i,j)*(WW(i,j+1)-WW(i,j-1)) )*dt/(2*dx);
    end
end



while 1
    

  
for i=2:num_x-1
    for j=2:num_y-1
       
        S(i,j)=0.25*(SS(i+1,j)+SS(i-1,j)+SS(i,j+1)+SS(i,j-1))-dx*dy*WW(i,j);
        
       dss(i,j)=abs(S(i,j)-SS(i,j));
        
    end
end

 maxa=dss(1,1);
 

		for i=2:num_x-1

			for j=2:num_y-1

				if dss(i,j)>maxa
					maxa=dss(i,j);
                end

			   
            end
        end




		




    if maxa<e
        break;
    end


SS=S;  
    
end

end


for i=2:num_x
    for j=2:num_y-1
        u(i,j)=(S(i,j+1)-S(i,j-1))/(2*dy);
        v(i,j)=(-S(i,j)+S(i-1,j))/(2*dx);
    end
end



quiver(y,y,u(90,:),v(90,:));
xlabel('u');
ylabel('y');



