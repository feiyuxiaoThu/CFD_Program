function x=GS(A,b,xi,eps,N)
%xΪ������Ľ�AΪϵ������bΪ������x0Ϊ������ֵepsΪ���N���޶��ĵ�������
%����Ҫ��A�ֽ�Ϊ�������Ǿ���
L=triu(A)-A;
U=tril(A)-A;
D=A+L+U;
Bs=U/(D-L);
fs=b/(D-L);
%�õ�������ʽBsΪ������fsΪ������
i=0;
con=0;
%����con��������¼�������Ƿ�����
while i<N
    i=i+1;
    x=Bs*xi+fs;
    for j=1:length(b)
        il(i,j)=x(j);
    end
    if norm(x-xi)<eps
        con=1;
        break
    end
    xi=x;
end
 
 