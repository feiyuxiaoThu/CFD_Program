function [Hx,Hy] = Velocity(u,v,dxu,dxs,dyv,dys,Re)
n = length(u);
% X-Momentum
ue2 = (1/2.*(u(4:n,2:n-1)+u(3:n-1,2:n-1))).^2;
uw2 = (1/2.*(u(2:n-2,2:n-1)+u(3:n-1,2:n-1))).^2;
du2dx = (ue2-uw2)./dxs;
un = 1/2.*(u(3:n-1,2:n-1) + u(3:n-1,3:n));
vn = 1/2.*(v(2:n-2,3:n) + v(3:n-1,3:n));
us = 1/2.*(u(3:n-1,1:n-2) + u(3:n-1,2:n-1));
vs = 1/2.*(v(2:n-2,2:n-1) + v(3:n-1,2:n-1));
duvdy = (un.*vn - us.*vs)./dyv;
d2udx2 = (u(4:n,2:n-1) - 2.*u(3:n-1,2:n-1) +u(2:n-2,2:n-1))./(dxu*dxs);
d2udy2 = (u(3:n-1,3:n) - 2.*u(3:n-1,2:n-1) +u(3:n-1,1:n-2))./(dys*dyv);
Hx(3:n-1,2:n-1) = -du2dx - duvdy + (1/Re).*(d2udx2+ d2udy2);
% Y-Momentum
vn2 = (1/2.*(v(2:n-1,4:n) + v(2:n-1,3:n-1))).^2;
vs2 = (1/2.*(v(2:n-1,3:n-1) + v(2:n-1,2:n-2))).^2;
dv2dy = (vn2 - vs2)./dys;
ue = 1/2.*(u(3:n,3:n-1) + u(3:n,2:n-2));
ve = 1/2.*(v(2:n-1,3:n-1) + v(3:n,3:n-1));
uw = 1/2.*(u(2:n-1,3:n-1) + u(2:n-1,2:n-2));
vw = 1/2.*(v(1:n-2,3:n-1) + v(2:n-1,3:n-1));
duvdx = (ue.*ve - uw.*vw)./dxu;
d2vdx2 = (v(3:n,3:n-1) - 2.*v(2:n-1,3:n-1) +v(1:n-2,3:n-1))./(dxs*dxu);
d2vdy2 = (v(2:n-1,4:n) - 2.*v(2:n-1,3:n-1) +v(2:n-1,2:n-2))./(dyv*dys);
Hy(2:n-1,3:n-1) = -dv2dy - duvdx + (1/Re).*(d2vdx2+ d2vdy2);