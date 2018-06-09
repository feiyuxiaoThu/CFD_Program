%final
clear;
D=1;
L=2;
U=1;
t_total=5;

Re=200;

dx=0.05;
dy=dx;
dt=0.005;
num_x=L/dx+1;
num_y=D/dy+1;
num_t=t_total/dt+1;

psi=zeros(num_x,num_y);
ppsi=zeros(num_x,num_y);
w=zeros(num_x,num_y);
ww=zeros(num_x,num_y);
u=zeros(num_x,num_y);
v=zeros(num_x,num_y);

times=40;

for i=1:num_x
    for j=2:num_y-1
        u(i,j)=U;
    end
end


%壁面边界条件

for i=1:num_x
    psi(i,1)=0;
    psi(i,num_y)=U*D;
end

% for i=1:num_x
%     w(i,1)=2*psi(i,2)/(dy*dy);
%     w(i,num_y)=(2*psi(i,num_y-1)-2*U*D)/(dy*dy);  %改变
% end

%进口
for i=1:num_x
    for j=1:num_y
    psi(i,j)=U*(j-1)*dy;
    end
end

for j=2:num_y-1
    w(1,j)=2*(psi(2,j)-psi(1,j))/(dx*dx)+(psi(1,j+1)-2*psi(1,j)+psi(1,j-1))/(dy*dy);  %改变
end

%出口
for j=2:num_y-1
    psi(num_x,j)=2*psi(num_x-1,j)-psi(num_x-2,j);    %改变
    w(num_x,j)=w(num_x-1,j);   %改变
end

for t=0:dt:t_total
   
    for i=2:num_x-1
        for j=2:num_y-1
            ww(i,j)=w(i,j)-dt*u(i,j)*(w(i+1,j)-w(i-1,j))/(2*dx)-dt*v(i,j)*(w(i,j+1)-w(i,j-1))/(2*dy)+(dt/Re)*((w(i+1,j)-2*w(i,j)+w(i-1,j))/(dx*dx)+(w(i,j+1)-2*w(i,j)+w(i,j-1))/(dy*dy));
        end
    end
     for i=2:num_x-1
        for j=2:num_y-1
            w(i,j)=ww(i,j);
        end
     end
     
    ppsi=psi;
    for s=1:times
        
        for i=2:num_x-1
            for j=2:num_y-1
                ppsi(i,j)=0.25*(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1))-dx*dy*w(i,j);
            end
        end
        
        for i=2:num_x-1
            for j=2:num_y-1
            psi(i,j)=ppsi(i,j);
            end
        end
        
    end
    
    
    %速度场
    for i=2:num_x-1
        for j=2:num_y-1
           u(i,j)=(psi(i,j+1)-psi(i,j-1))/(2*dy);
           v(i,j)=-(psi(i+1,j)-psi(i-1,j))/(2*dx);
        end
    end
     
    

for i=1:num_x
    w(i,1)=2*psi(i,2)/(dy*dy);
    w(i,num_y)=(2*psi(i,num_y-1)-2*U*D)/(dy*dy);  %改变
end

%进口


for j=2:num_y-1
    w(1,j)=2*(psi(2,j)-psi(1,j))/(dx*dx)+(psi(1,j+1)-2*psi(1,j)+psi(1,j-1))/(dy*dy);  %改变
end

%出口
for j=2:num_y-1
    psi(num_x,j)=2*psi(num_x-1,j)-psi(num_x-2,j);    %改变
    w(num_x,j)=w(num_x-1,j);   %改变
end
    
end

d=zeros(num_x,num_y);
e=0.01;
while 1
    for i=2:num_x-1
            for j=2:num_y-1
                ppsi(i,j)=0.25*(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1))-dx*dy*w(i,j);
                d(i,j)=abs(ppsi(i,j)-psi(i,j));
            end
    end
    
    for i=2:num_x-1
            for j=2:num_y-1
            psi(i,j)=ppsi(i,j);
            end
        end
        
        
    
    maxa=d(1,1);
 

		for i=2:num_x-1

			for j=2:num_y-1

				if d(i,j)>maxa
					maxa=d(i,j);
                end

			   
            end
        end




		

 


    if maxa<e
        break;
    end
    
end
    
plot(u(10,:));
xlabel('y');
ylabel('u');
mesh(u);
xlabel('y');
ylabel('x');
zlabel('u');