function [x]=chase(a,c,d,b)%����a�洢���Ǿ���A�����Խ���Ԫ�أ�c��d�洢���Խ����ϱ��±ߴ���Ϊ1��Ԫ��
 
    n=length(a);
    n1=length(c);
    n2=length(d);
    %������
    if n1~=n2%�洢���������ά������
        error('MATLAB:Crout:�������ԽǾ��󣬲���������Ԫ�ظ�������.');
    elseif n~=n1+1
        error('MATLAB:Crout:�������ԽǾ��󣬲���������Ԫ�ظ�������.');
    end
   
    %��ʼ��
   
    p=1:n;
    q=1:n-1;
    x=1:n;
    y=1:n;
   
    %׷�Ϸ���������
    p(1)=a(1);
    for i=1:n-1
        q(i)=c(i)/p(i);
        p(i+1)=a(i+1)-d(i)*q(i);%d���±��Ϊ1��n-1
    end
    %����y
    y(1)=b(1)/p(1);%��x�洢y
    for i=2:n
        y(i)=(b(i)-d(i-1)*y(i-1))/p(i);
    end
    %����x
    x(n)=y(n);
    for i=(n-1):-1:1
        x(i)=y(i)-q(i)*x(i+1);
    end
   
end %end of function