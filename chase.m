function [x]=chase(a,c,d,b)%数组a存储三角矩阵A的主对角线元素，c、d存储主对角线上边下边带宽为1的元素
 
    n=length(a);
    n1=length(c);
    n2=length(d);
    %错误检查
    if n1~=n2%存储矩阵的数组维数错误
        error('MATLAB:Crout:不是三对角矩阵，参数数组中元素个数错误.');
    elseif n~=n1+1
        error('MATLAB:Crout:不是三对角矩阵，参数数组中元素个数错误.');
    end
   
    %初始化
   
    p=1:n;
    q=1:n-1;
    x=1:n;
    y=1:n;
   
    %追赶法程序主体
    p(1)=a(1);
    for i=1:n-1
        q(i)=c(i)/p(i);
        p(i+1)=a(i+1)-d(i)*q(i);%d的下标改为1到n-1
    end
    %正解y
    y(1)=b(1)/p(1);%用x存储y
    for i=2:n
        y(i)=(b(i)-d(i-1)*y(i-1))/p(i);
    end
    %倒解x
    x(n)=y(n);
    for i=(n-1):-1:1
        x(i)=y(i)-q(i)*x(i+1);
    end
   
end %end of function