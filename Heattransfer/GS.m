function x=GS(A,b,xi,eps,N)
%x为方程组的解A为系数矩阵b为常数项x0为迭代初值eps为误差N是限定的迭代次数
%首先要将A分解为上下三角矩阵
L=triu(A)-A;
U=tril(A)-A;
D=A+L+U;
Bs=U/(D-L);
fs=b/(D-L);
%得到迭代格式Bs为迭代阵fs为常向量
i=0;
con=0;
%其中con是用来记录计算结果是否收敛
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
 
 