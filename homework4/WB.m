%������Բ�������
%Warming-Beam��ʽ
clear all;
Mx=100;
CFL=1.5;
cc=CFL*CFL;
a=1;
t_total=1;
dt_record=0.005;

Dx=0.5-(-0.5);
dx=Dx/Mx;

dt=dx*CFL/a;

Mt=ceil(t_total/dt);  %ʱ������

Ux=zeros(1,Mx);
Uxx=zeros(1,Mx);
y=zeros(1,Mx);

x=zeros(1,Mx);

for i=1:1:Mx
    x(i)=-0.5+dx*(i-1);
end




%��ʼ��
for i=1:1:Mx
    Ux(1,i)=fx(x(i));
    y(1,i)=fx(x(i)-rem(t_total,Dx));
end

Uxx(1,:)=Ux(1,:);

for i=1:1:Mt
    
    u0=Uxx(1,Mx-2);
    um1=Uxx(1,Mx-1);
    Ux(1,1)=(1-1.5*CFL+cc/2)*Uxx(1,1)+(2*CFL-cc)*u0+(cc/2-CFL/2)*um1;
    Ux(1,2)=(1-1.5*CFL+cc/2)*Uxx(1,2)+(2*CFL-cc)*Uxx(1,1)+(cc/2-CFL/2)*u0;
    for ii=3:1:Mx;
        
        Ux(1,ii)=(1-1.5*CFL+cc/2)*Uxx(1,ii)+(2*CFL-cc)*Uxx(1,ii-1)+(cc/2-CFL/2)*Uxx(1,ii-2);
    
    end
    
    
    
    
    Uxx(1,:)=Ux(1,:);
    
%     plot(Uxx(1,:));
%     hold on;
    
end

    
 plot(y);
 hold on;
 plot(Uxx);
     
   %axis([-0.5 0.5 0 t_total 0 1]);
   title('Warming-Beam��ʽ,t_total=10');
   xlabel('x');
   ylabel('t');


