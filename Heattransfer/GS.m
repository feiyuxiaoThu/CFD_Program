function  x=GS(a,b)

e=1e-4;

n=length(b);

N=10000;

x=zeros(n,1);

tt=zeros(n,1);

for kk=1:N

    sum=0;

    tt(1:n)=x(1:n);

    for i=1:n      

    x(i)=(b(i)-a(i,1:(i-1))*x(1:(i-1))-a(i,(i+1):n)*tt((i+1):n))/a(i,i);

    end

    for i=1:n

        sum=sum+(tt(i)-x(i))^2;

    end

    if sqrt(sum)<e  

        break;

    end

end

if kk==N warning('未能找到近似解');

end

 