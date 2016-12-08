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
        sum = sum + A(i,j)*x(j);
    end
    x(i) = (b(i) - sum)/A(i,i);
end
x = x';


end
