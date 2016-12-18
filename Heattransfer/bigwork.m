% main.m
%拼装矩阵
clear;
L=0.1; %矩形边长
k=5; %热导率
tf1=80;
h1=50; %左边界的对流参数
tf2=30;
h2=100; %右边界的对流参数

N=30; %划分格子数目
d=L/N;

num=(N+1)*(N+1);

T=zeros(1,num);
A=zeros(num,num);
b=zeros(1,num);


%内部点
for i=2:N
    for j=2:N
        A((i-1)*(N+1)+j,(i-1)*(N+1)+j+1)=0.25;
         A((i-1)*(N+1)+j,(i-1)*(N+1)+j-1)=0.25;
          A((i-1)*(N+1)+j,(i-2)*(N+1)+j)=0.25;
           A((i-1)*(N+1)+j,i*(N+1)+j)=0.25;
    end
end

%上边界绝热
for i=2:N
    A((i-1)*(N+1)+N+1,(i-2)*(N+1)+N+1)=0.25;
    A((i-1)*(N+1)+N+1,i*(N+1)+N+1)=0.25;
    A((i-1)*(N+1)+N+1,(i-1)*(N+1)+N)=0.5;
end

%下边界恒定温度
for i=2:N
    b((i-1)*(N+1)+1)=120;
end

%左边界对流
for i=2:N
    A(i,N+1+i)=2/(h2*d/k +4);
    A(i,i+1)=1/(h2*d/k +4);
    A(i,i-1)=1/(h2*d/k +4);
    b(i)=(h2*d*tf2/k)/(h2*d/k +4);
end

%右边界对流

for i=2:N
    A(N*(N+1)+i,(N-1)*(N+1)+i)=2/(h1*d/k +4);
    A(N*(N+1)+i,N*(N+1)+i+1)=1/(h1*d/k +4);
    A(N*(N+1)+i,N*(N+1)+i-1)=1/(h1*d/k +4);
    b(N*(N+1)+i)=(h1*d*tf1/k)/(h1*d/k +4);
end

%角点
%t11
A(1,(N+1)+1)=1/(h2*d/k + 2);
A(1,2)=1/(h2*d/k +2 );
b(1)=(d*h2*tf2/k)/(h2*d/k +2);
%t1(N+1)
A(N+1,N)=1/(h2*d/k + 2);
A(N+1,(N+1)+N+1)=1/(h2*d/k + 2);
b(N+1)=(d*h2*tf2/k)/(h2*d/k +2);
%t(N+1)1
A(N*(N+1)+1,(N-1)*(N+1)+1)=1/(h1*d/k + 2);
A(N*(N+1)+1,N*(N+1)+2)=1/(h1*d/k + 2);
b(N*(N+1)+1)=(d*h1*tf1/k)/(h1*d/k +2);
%t(N+1)(N+1)
A(N*(N+1)+N+1,(N-1)*(N+1)+N+1)=1/(h1*d/k + 2);
A(N*(N+1)+N+1,N*(N+1)+N)=1/(h1*d/k + 2);
b(N*(N+1)+N+1)=(d*h1*tf1/k)/(h1*d/k +2);



 AA=eye(num)-A;

 t0=cputime;
% x=Jocabi(AA,b);
 %x = Gauss(AA,b);
 x=GS(AA,b);
 %x = SOR( AA, b, T, 1.1, 1e-4, 5000 );
 t1=cputime-t0;

t=zeros(N+1,N+1);
for i=1:N+1
    for j=1:N+1
        t(i,j)=x((i-1)*(N+1)+j);
    end
end

figure ;
mesh(t');
xlabel('y');
ylabel('x');
zlabel('T');
title('温度场分布');

%绘制图形
p_l=zeros(1,N-1);
p_r=zeros(1,N-1);
p_d=zeros(1,N-1);
for i=1:N-1
    p_l(i)=h2*(tf2-t(1,i+1))*d;
    p_r(i)=h1*(tf1-t(N+1,i+1))*d;
    p_d(i)=k*(t(i+1,1)-t(i+1,2));
end

%
figure;
 plot(p_l);
 hold on;
 plot(p_r);
 hold on;
 plot(p_d);
 legend('左边界热流','右边界热流','下边界热流');
 title('边界热流分布');


%Guass消元法
function x = Gauss(A,b)
%Gauss Gauss消元法解方程
[~,n] = size(A);
for i = 1:(n-1);
  for j = (i+1):n;
      l = A(j,i)/A(i,i);
      A(j,:) = A(j,:) - l*A(i,:);
      b(j)   = b(j) - l*b(i);
  end
end
x(1:n) = 0;
x(n) = b(n)/A(n,n);
for i = (n-1):-1:1;
    sum = 0;
    for j = (i+1):n;
        sum = sum + A(i,j)*x(
        j);
    end
    x(i) = (b(i) - sum)/A(i,i);
end
x = x';



end

%高斯迭代法
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

%G-S迭代法
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


 %SOR超松弛迭代法

function x = SOR( A, b, xold, omega, TOL, Nmax )

%SOR          approximate the solution of the linear system Ax = b by applying
%             the SOR method (successive over-relaxation)
%
%     calling sequences:
%             x = sor ( A, b, xold, omega, TOL, Nmax )
%             sor ( A, b, xold, omega, TOL, Nmax )
%
%     inputs:
%             A       coefficient matrix for linear system - must be a
%                     square matrix
%             b       right-hand side vector for linear system
%             xold    vector containing initial guess for solution of
%                     linear system
%             omega   relaxation parameter
%             TOL     convergence tolerance - applied to maximum norm of
%                     difference between successive approximations
%             NMax    maximum number of iterations to be performed
%
%     output:
%             x       approximate solution of linear system
%
%     NOTE:
%             if SOR is called with no output arguments, the
%             iteration number and the current approximation to the
%             solution are displayed
%
%             if the maximum number of iterations is exceeded, a meesage
%             to this effect will be displayed and the current approximation
%             will be returned in the output value
%

n = length ( b );
[r c] = size ( A );
if ( c ~= n )
   disp ( 'sor error: matrix dimensions and vector dimension not compatible' )
   return
end;
xnew = zeros ( 1, n );

if ( nargout == 0 )
   s = sprintf ( '%3d \t %10f ', 0, xold(1) );
   for j = 2 : n
	   s = sprintf ( '%s%10f ', s, xold(j) );
   end;
   disp ( s );
end;

for its = 1 : Nmax
    xnew(1) = ( 1 - omega ) * xold(1) + omega * ( b(1) - sum ( A(1,2:n) .* xold(2:n) ) ) / A(1,1);
	for i = 2 : n-1
	    xnew(i) = ( 1 - omega ) * xold(i) + omega * ( b(i) - sum ( A(i,1:i-1) .* xnew(1:i-1) ) - sum ( A(i,i+1:n) .* xold(i+1:n) ) ) / A(i,i);
	end;
	xnew(n) = ( 1 - omega ) * xold(n) + omega * ( b(n) - sum ( A(n,1:n-1) .* xnew(1:n-1) ) ) / A(n,n);

    if ( nargout == 0 )
	   s = sprintf ( '%3d \t %10f ', its, xnew(1) );
	   for j = 2 : n
	       s = sprintf ( '%s%10f ', s, xnew(j) );
	   end;
	   disp ( s );
	end;

    conv = max ( abs ( xnew - xold ) );
	if ( conv < TOL )
	   x = xnew;
	   return
	else
	   xold = xnew;
	end;
end;
disp ( 'sor error: maximum number of iterations exceeded' );
if ( nargout == 1 ) x = xnew; end;
