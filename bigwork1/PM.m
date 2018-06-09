clear;
A = importdata('Data.txt');
n0=size(A)/2;
N=n0(1);

n=N-2;

x=zeros(1,n+2);
y=zeros(1,n+2);



for i=1:n+2
    x(i)=A(i);
    y(i)=A(N+i);
end
pi=3.1415926;
%来流
U=5;
Angle_a=5*pi/180; %十度仰角
%n+1个面元
xm=zeros(1,n+1);
ym=zeros(1,n+1);
An=zeros(n+1,n+1);
Bn=zeros(n+1,n+1);
At=zeros(n+1,n+1);
Bt=zeros(n+1,n+1);
r=zeros(n+1,n+1);
rr=zeros(n+1,n+1);
Vn=zeros(1,n+1);
Vt=zeros(1,n+1);

thlta=zeros(1,n+1);
belta=zeros(n+1,n+1);

for i=1:n+1
    xm(i)=0.5*(x(i)+x(i+1));
    ym(i)=0.5*(y(i)+y(i+1));
    thlta(i)=atan((y(i+1)-y(i))/(x(i+1)-x(i)));
end

for i=1:n+1
    for j=1:n+1
    r(i,j)=((xm(i)-x(j))^2+(ym(i)-y(j))^2)^0.5;
    rr(i,j)=((xm(i)-x(j+1))^2+(ym(i)-y(j+1))^2)^0.5;
    belta(i,j)=atan((ym(i)-y(j+1))/(xm(i)-x(j+1)))-atan((ym(i)-y(j))/(xm(i)-x(j)));
    end
end

for i=1:n+1
    for j=1:n+1
        if i==j
            An(i,j)=0.5;
           At(i,j)=0;
        else   
        An(i,j)=(0.5/pi)*(sin(thlta(i)-thlta(j))*log(rr(i,j)/r(i,j))+cos(thlta(i)-thlta(j))*belta(i,j));
        At(i,j)=(0.5/pi)*(sin(thlta(i)-thlta(j))*belta(i,j)-cos(thlta(i)-thlta(j))*log(rr(i,j)/r(i,j)));
        end
        
        Bn(i,j)=-At(i,j);
        Bt(i,j)=An(i,j);
    end
end

AA=zeros(n+2:n+2);
for i=1:n+1
    for j=1:n+1
        AA(i,j)=An(i,j);
    end
end

for i=1:n+1
    for j=1:n+1
    AA(i,n+2)=AA(i,n+2)+Bn(i,j);
    end
end

for i=1:n+1
    AA(n+2,i)=At(1,j)+At(n+1,j);
    AA(n+2,n+2)=AA(n+2,n+2)+Bt(1,j)+Bt(n+1,j);
end


%b
b=zeros(1,n+2);
for i=1:n+1
    b(i)=-U*sin(Angle_a-thlta(i));
end

b(n+2)=-U*cos(Angle_a-thlta(1))-U*cos(Angle_a-thlta(n+1));


q=zeros(1,n+2);

q=AA\b';

for i=1:n+1
    for j=1:n+1
    Vt(i)=Vt(i)+At(i,j)*q(j)+q(n+2)*Bt(i,j);
    end
end

for i=1:n+1
    Vt(i)=Vt(i)+U*cos(Angle_a-thlta(i));
end

%求解
Cv=zeros(1,n+1);
Cp=zeros(1,n+1);
for i=1:n+1
    Cv(i)=Vt(i)/U;
    Cp(i)=1-(Vt(i)/U)^2;
end

plot(Cp);
    