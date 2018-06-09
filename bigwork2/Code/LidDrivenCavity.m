clc
clear
%% 
%����
wt0 = cputime;
Re = 1000;  %Reynolds number
dt = 0.0005;  %ʱ�䲽��
nc = 68;  %����������ĸ����
update = 100;
tmax = 20;  %������ʱ�䣨Poisson���̵���������
ipmax = 5000;  %������������Poisson���̵���������
U0 = 1;    %ƽ���ƶ����ٶ�
n = nc+2;    
L = 1;
nx = nc;
ny = nc;
d_cell = L/nc;
dxu = d_cell;dxs = d_cell;dyv = d_cell;dys = d_cell;
x = 0:1/(nc-1):1;
y = 0:1/(nc-1):1;
%% 
%�ٶ���ʼ��
p = zeros(nc+2,nc+2);
u = zeros(nc+2,nc+2);
v = zeros(nc+2,nc+2);
us = u;
vs = v;
ps = p;
p_old = p;
u_old = u;
v_old = v;
S(1:n-1,1:n-1) = 0;
%% 
% Poisson���̵�ϵ��
As = zeros(n,n);
As(:,3:n-1) = dxu/dys;
An = zeros(n,n);
An(:,2:n-2) = dxu/dys;
Aw = zeros(n,n);
Aw(3:n-1,:) = dyv/dxs;
Ae = zeros(n,n);
Ae(2:n-2,:) = dyv/dxs;
Ap = -(Aw + Ae + As + An);
%����������
iter.p = 2;
iter.o = 1;
%�����в�
res.p(1) = 0;
tp = 0;
%%
%�������
for t = 0:dt:tmax
% �ü���ĳ�ʼ��ѹ���ֲ��õ��²���ٶȳ�
[u,v] = Boundaryconditions(u,v,U0);   %�߽�����
[Hx,Hy] = Velocity(u,v,dxu,dxs,dyv,dys,Re);  %��⶯������
if iter.o == 1
us(2:n-1,1:n-1) = u(2:n-1,1:n-1) +dt*Hx(2:n-1,1:n-1);
vs(1:n-1,2:n-1) = v(1:n-1,2:n-1) +dt*Hy(1:n-1,2:n-1);
else
[Hx_old,Hy_old] =Velocity(u_old,v_old,dxu,dxs,dyv,dys,Re);
us(2:n-1,1:n-1) = u(2:n-1,1:n-1) ...
+ dt/2*(3*Hx(2:n-1,1:n-1) -Hx_old(2:n-1,1:n-1));
vs(1:n-1,2:n-1) = v(1:n-1,2:n-1) ...
+ dt/2*(3*Hy(1:n-1,2:n-1) -Hy_old(1:n-1,2:n-1));
end
u_tp = us;
v_tp = vs;
% Poisson����
for i = 2:n-1
for j = 2:n-1
S(i,j) = 1/dt*((us(i+1,j) -us(i,j))*dyv + (vs(i,j+1) - vs(i,j))*dxu);  %�ұ�ϵ��
end
end
tpp = cputime;
for c = 1:ipmax
iter.p = iter.p + 1;
 %���ɷ������
[p,res.p(iter.p)] =PoissonGS(p,n,Ap,An,As,Ae,Aw,S);
if res.p(iter.p) < 1e-5   %ĳ�������ϵĵ������������������Ҷ���Ϊ����ѹ����ľ���ֵ�͵Ŀ���
break
end
end
tp = tp + (cputime - tpp);
% �����ٶȳ�
u(2:n-1,2:n-1) = us(2:n-1,2:n-1) ...
- dt/dxs*(p(2:n-1,2:n-1) - p(1:n-2,2:n-1));
v(2:n-1,2:n-1) = vs(2:n-1,2:n-1) ...
- dt/dys*(p(2:n-1,2:n-1) - p(2:n-1,1:n-2));
err.p(iter.o) = norm(p - p_old);
err.u(iter.o) = norm(u - u_old);
err.v(iter.o) = norm(v - v_old);
if err.u(iter.o) < 1E-6 && err.v(iter.o) < 1E-6
beep
wt = cputime - wt0;
break
end

p = p - p(5,5); % normalize pressure
p_old = p;
u_old = u;
v_old = v;

iter.o = iter.o + 1;
end

%%
%���ӻ�
uc = 0.5*(u(2:end-1,2:end-1) + u(3:end,2:end-1));
vc = 0.5*(v(2:end-1,2:end-1) + v(2:end-1,3:end));
[X,Y] = meshgrid(x,y);
figure
streamslice(X,Y,uc',vc',2)
axis([0,L,0,L])
title('Lid Driven Cavity ,Velocity');
xlabel('X')
ylabel('Y')

figure
quiver(uc',vc',0.6,'k-')
       axis equal
        title('Lid Driven Cavity,P ');
        hold on 
        pcolor(p(2:end-1,2:end-1)');
        colormap(jet)
        colorbar
        shading interp  
        drawnow;
        hold off

