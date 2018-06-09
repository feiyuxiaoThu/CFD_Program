%Aerofoil:
clear;
Aerofoil = 0012;

%chord:
c = 1;

%Thickness:
t = c*Aerofoil/100;

%Half the mesh resolution - 1:
Npanel = 90;

%This is to avoid any issues when one panel is taken out:
if Npanel < 8
    
    Npanel = 8;
    
end

%In angular terms,
delta_beta = pi/(Npanel);
beta = 0:delta_beta:2*pi;
%2*pi implies that top and bottom are being simulated.

%Corresponding x-coordinates:
x = c/2*(1 - cos(beta));
x_t = x(1:Npanel+1);
x_b = x(Npanel+1:2*Npanel);

%x_t is for the top side and goes from the leading edge to the trailing
%edge. x_b is for the bottom and goes from the trailing edge to the leading
%edge. This is because of the setup of the matrices defined at the bottom
%of page 282 in the Katz and Plotkin book which incorporate the Kutta
%condition into the simulation.

%Calculate the z-coordinates (top half only initially):
z_t = 5*t*(0.2969*sqrt(x_t/c) + (-0.1260)*(x_t/c) + (-0.3516)*(x_t/c).^2 + 0.2843*(x_t/c).^3 + (-0.1036)*(x_t/c).^4);
z_b = -(z_t(end:-1:2));

%Panel corner coordinates:
x = horzcat(x_b,x_t);
z = horzcat(z_b,z_t);
z(1) = 0;
z(end) = 0;



fid = fopen('Datax.txt','wt');%数据保存在你当前的文件夹下，文件名为Data.txt
Temp1 = x;
fprintf(fid,'%f\n',Temp1);
fclose(fid);

fid = fopen('Datay.txt','wt');%数据保存在你当前的文件夹下，文件名为Data.txt
Temp2 = z;
fprintf(fid,'%f\n',Temp2);
fclose(fid);


%Show the aerofoil:
figure
plot(x,z);