function y = myfun(x,i) 

b(1)=a1(i);
b(2)=a2(i);
b(3)=a3(i);
b(4)=a4(i);
b(5)=a5(i);
c=b\A;
y = (c(1)*sin(-3*x)+c(2)*sin(-2*x)+c(3)*sin(-x)+c(4)*sin(x)+a6(i)*sin(2*x)-x).^2;
