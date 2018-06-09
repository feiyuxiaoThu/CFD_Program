clear;
A=[1,1,1,1,1;-3,-2,-1,0,1;9,4,1,0,1;-27,-8,-1,0,1;81,16,1,0,1];
dx=0.0001;
a=-0.3:dx:0.3;
N=6001;
b=zeros(1,5);
y=zeros(1,N);
for i=1:N
    b(1)=-a(i);
    b(2)=1-2*a(i);
    b(3)=-4*a(i);
    b(4)=-8*a(i);
    b(5)=-16*a(i);
    c=inv(A)*b';
    
    
     F = @(x)exp(9*(3.1416-x)).*(c(1)*sin(-3*x)+c(2)*sin(-2*x)+c(3)*sin(-x)+c(5)*sin(x)+a(i)*sin(2*x)-x).^2;
     y(i) = quad(F,0,3.1416); 
    
    
end

plot(a,y);
xlabel('a');
ylabel('E');