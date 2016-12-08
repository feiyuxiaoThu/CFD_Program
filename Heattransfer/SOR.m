function x = SOR(A,b,w,m)
%sor sor法解方程
%w为松弛引子,m为迭代次数
 [~,n] = size(A);
 D = diag(diag(A));
 L = D - tril(A);
 U = D - triu(A);
 BW= (D - w*L)\((1-w)*D+w*U);
 FW= w*((D - w*L)\b);
 x(1:n) = 0;
 x = x';
 for i = 1:m;
     x = BW*x + FW;
 end

end
