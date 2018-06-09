function x = Tri2(a,b,c,f)
n = length(b);
beta(1) = b(1);
gamma(1) = f(1)/beta(1);
for q = 2:n
beta(q) = b(q) - a(q)*c(q-1)/beta(q-1);
gamma(q) = (-a(q)*gamma(q-1) +f(q))/beta(q);
end
x(n) = gamma(n);
for q = n-1:-1:1
x(q) = gamma(q) - x(q+1)*c(q)/beta(q);
end