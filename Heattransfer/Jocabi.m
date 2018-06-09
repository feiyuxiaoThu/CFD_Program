function  x=Jocabi(a,b)

e=1e-4;

n=length(b);

N=5000;

x=zeros(n,1);

y=zeros(n,1);

for k=1:N

    sum=0;

    for i=1:n

    y(i)=(b(i)-a(i,1:n)*x(1:n)+a(i,i)*x(i))/a(i,i);

    end

    for i=1:n

        sum=sum+(y(i)-x(i))^2;

    end

    if sqrt(sum)<e

        k      

        break;    

    else

        for i=1:n

        x(i)=y(i);

        end

    end

end

if k==N warning('未能找到近似解');

end