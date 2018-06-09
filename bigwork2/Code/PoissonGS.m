function [phi,res] =PoissonGS(phi,n,Ap,An,As,Ae,Aw,S)
%%Gauss-Sedialµü´ú·¨
phi_new=phi;
for i = 2:n-1
    for j=2:n-1
    phi(i,j)=(1/Ap(i,j))*(S(i,j)-As(i,j)*phi(i,j-1)-Aw(i,j)*phi(i-1,j)-Ae(i,j)*phi(i+1,j)-An(i,j)*phi(i,j+1));
    end
end
res=sqrt(sum(sum(abs(phi_new-phi))));
