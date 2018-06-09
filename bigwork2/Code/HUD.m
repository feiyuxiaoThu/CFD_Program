function [] = HUD(t,wt0,iter,err,tstring)
io = 1:iter.o;
semilogy(io,err.p,io,err.u,io,err.v)
title(['Lid Driven Cavity Error, ',tstring])
legend('pressure','u-velocity','v-velocity')
xlabel('Iteration')
ylabel('Error')
drawnow expose
wt = cputime - wt0;
clc
fprintf('Solution time (s): %4.3e\n',t);
fprintf('Wall time (s) : %4.3e\n',wt);
fprintf('Speed : %4.3e\n\n',t/wt);
fprintf('Errors\n')
fprintf(' p:%4.3e\n',err.p(iter.o));
fprintf(' u:%4.3e\n',err.u(iter.o));
fprintf(' v:%4.3e\n',err.v(iter.o));