function [u,v] = Boundaryconditions(u,v,U0)
% West
v(1,:) = -v(2,:);
u(2,:) = 0;
% East
v(end,:) = -v(end-1,:);
u(end,:) = 0;
% South
u(:,1) = - u(:,2);
v(:,2) = 0;
% North
u(:,end) = (2*U0 - u(:,end-1));
v(:,end) = 0;