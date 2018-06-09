%bigwork1
clear ;
%读入翼型表面数据
fid=fopen('xyData.dat');
B=textscan(fid,'%f %f');%把每一列的数据读入到读入到单元数组B中
C=[B{1} B{2}];          %从单元数组B中提取每列数据赋值给矩阵C
n=max(size(C));         %确定读入数据的坐标数目
x0=C(:,1);y0=C(:,2);      %赋值
Dy=15;
dy=0.01;
dyy=0.01;
pi=3.1415926;

dn=0.1;
de=0.1;

dr=0.01;

thlta=zeros(1,n);

numrx=60;
num_x=2*numrx+n;
num_y=60;

x=zeros(num_x,num_y);
y=zeros(num_x,num_y);

x(numrx,1)=1.008;
y(numrx,1)=0;
x(numrx+n+1,1)=1.008;
y(numrx+n+1,1)=0;

%翼型表面
for i=numrx+1:numrx+n
    x(i,1)=x0(i-numrx);
    y(i,1)=y0(i-numrx);
end

%尾流区的边界
for i=numrx-1:-1:1
    x(i,1)=x(i+1,1)+(numrx-i)*dy;
    y(i,1)=0;
end
for i=numrx+n+2:1:num_x
    x(i,1)=x(num_x+1-i,1);
    y(i,1)=0;
end

x(numrx,num_y)=x(numrx,1);
y(numrx,num_y)=Dy;
x(numrx+n+1,num_y)=x(numrx,1);
y(numrx+1,num_y)=-Dy;
%尾流区内部
for j=2:num_y
    for i=numrx:-1:1
        x(i,j)=x(i+1,1)+(numrx-i)*dy;
        y(i,j)=dy*(j-1)+y(i,j-1);
    end
end

for j=2:num_y
    for i=numrx+n:1:num_x
        x(i,j)=x(num_x+1-i,j);
        y(i,j)=-dy*(j-1)+y(i,j-1);
    end
end


for i=1:n
    if(i<41)
        thlta(i)=i*pi/(2*41);
    else
        thlta(i)=(n-i)*pi/(2*41);
    end
end


%左边半圆形区域区域
for i=numrx+1:numrx+40
    
        %x(i,num_y)=x(numrx,1)-Dy*cos(thlta(i-numrx));
        x(i,num_y)=x(numrx,num_y)-0.375*(i-numrx);
         y(i,num_y)=sqrt(Dy*Dy-(x(i,num_y)-1)*(x(i,num_y)-1));
    
end
 x(numrx+41,num_y)=-14;
 y(numrx+41,num_y)=0;
for i=numrx+42:numrx+n
    
    
    
      % x(i,num_y)=x(numrx+n,1)-Dy*cos(thlta(i-numrx));
      x(i,num_y)=x(n+numrx-i+numrx,num_y);
     y(i,num_y)=-sqrt(Dy*Dy-(x(i,num_y)-1)*(x(i,num_y)-1));
    
    
end


    
%圆形区域内部
for j=2:num_y
    for i=numrx:1:numrx+n
        
        if i<numrx+41
            x(i,j)=x(i,1)-(y(numrx,j)-y(numrx,1))*sin(thlta(i-numrx+1));
            y(i,j)=y(i,1)+(y(numrx,j)-y(numrx,1))*cos(thlta(i-numrx+1));
        else
            x(i,j)=x(i,1)-(y(numrx,j)-y(numrx,1))*sin(thlta(i-numrx));
            y(i,j)=y(i,1)-(y(numrx,j)-y(numrx,1))*cos(thlta(i-numrx));
        end
        
        y(numrx+40,j)=0;
    
    end       
            
    end


%生成网格
a=zeros(num_x,num_y);
b=zeros(num_x,num_y);
dx=zeros(num_x,num_y);
dy=zeros(num_x,num_y);

alph=zeros(num_x,num_y);
beta=zeros(num_x,num_y);
gama=zeros(num_x,num_y);
bw=zeros(num_x,num_y);
be=zeros(num_x,num_y);
bs=zeros(num_x,num_y);
bn=zeros(num_x,num_y);
bp=zeros(num_x,num_y);

cpx=zeros(num_x,num_y);
cpy=zeros(num_x,num_y);

s=0;
e=0.001;



while 1
    for i=1:num_x
        for j=1:num_y
             a(i,j)=x(i,j);
            b(i,j)=y(i,j);
        end
    end
    
    for i=2:num_x-1
        for j=2:num_y-1
          
            alph(i,j)=((x(i,j+1)-x(i,j-1))/(2*dn))^2+((y(i,j+1)-y(i,j-1))/(2*dn))^2;
            beta(i,j)=((x(i+1,j)-x(i-1,j)))*(x(i,j+1)-x(i,j-1))/(4*dn*de)+((y(i+1,j)-y(i-1,j)))*(y(i,j+1)-y(i,j-1))/(4*dn*de);
            gama(i,j)=( (x(i+1,j)-x(i-1,j))/(2*dn))^2+((y(i+1,j)-y(i-1,j))/(2*de))^2;

            bw(i,j)=alph(i,j)/(de*de);
            be(i,j)=bw(i,j);

            bs(i,j)=gama(i,j)/(dn*dn);
            bn(i,j)=bs(i,j);

            bp(i,j)=bw(i,j)+be(i,j)+bs(i,j)+bn(i,j);

            cpx(i,j)=-beta(i,j)*(x(i+1,j+1)-x(i+1,j-1)-x(i-1,j+1)+x(i-1,j-1))/(2*de*dn);
            cpy(i,j)=-beta(i,j)*(y(i+1,j+1)-y(i+1,j-1)-y(i-1,j+1)+y(i-1,j-1))/(2*de*dn);
        end
    end


    for i=2:num_x-1
        for j=2:num_y-1
            x(i,j)=(bw(i,j)*a(i-1,j)+be(i,j)*a(i+1,j)+bs(i,j)*a(i,j-1)+bn(i,j)*a(i,j+1)+cpx(i,j))/bp(i,j);
            y(i,j)=(bw(i,j)*b(i-1,j)+be(i,j)*b(i+1,j)+bs(i,j)*b(i,j-1)+bn(i,j)*b(i,j+1)+cpy(i,j))/bp(i,j);

            dx(i,j)=abs(a(i,j)-x(i,j));
            dy(i,j)=abs(b(i,j)-y(i,j));
        end
    end

    maxa=dx(1,1);
    maxb=dy(1,1);

		for i=2:num_x-1

			for j=2:num_y-1

				if dx(i,j)>=maxa

					maxa=dx(i,j);
                end

			    if dy(i,j)>=maxb

					maxb=dy(i,j);
                end
            end
        end




		s=s+1;




    if maxa<e && maxb<e
        break;
    end
end

  z=zeros(num_x,num_y);
 mesh(z,x,y);
%plot(x(:,num_y),y(:,num_y))
fid0=fopen('mu.txt');
BB=textscan(fid0,'%f %f %f');%把每一列的数据读入到读入到单元数组B中
%CCC=[BB%nn=max(size(CCC));  
%mu=CCC(:,1);
%x_c=CCC(:,2);z_c=CCC(:,3); 
mu=BB{1,1};
x_c=BB{1,2};
z_c=BB{1,3};
u_infty=1;
Uxw=zeros(num_x,num_y);
Uyw=zeros(num_x,num_y);

for i=2:num_x-1
        for j=2:num_y-1
    Uxw(i,j)=u_infty*cos(10*180/pi);
     Uyw(i,j)=u_infty*sin(10*180/pi);
        end
end


for k=2:num_x-1
        for j=2:num_y-1


    for i=1:179
        Uxw(k,j)=Uxw(k,j)+mu(i)*(x(k,j)-x_c(i))/((2*pi)*((x(k,j)-x_c(i))^2+(y(k,j)-z_c(i))^2))+mu(180)*(y(k,j)-z_c(i))/((2*pi)*((x(k,j)-x_c(i))^2+(y(k,j)-z_c(i))^2));
        Uyw(k,j)=Uyw(k,j)+mu(i)*(y(k,j)-z_c(i))/((2*pi)*((x(k,j)-x_c(i))^2+(y(k,j)-z_c(i))^2))-mu(180)*(x(k,j)-x_c(i))/((2*pi)*((x(k,j)-x_c(i))^2+(y(k,j)-z_c(i))^2));
        
    end
        end
        
        
end


 
fid = fopen('xy.txt', 'wt');
fprintf(fid, '%f %f\n', x,y);
 fclose(fid);
