% plot 学习
[X,Y]=meshgrid(-2:0.5:2,-2:0.5:2);  %生成坐标轴
Z=2*X.*exp(-X.^2-Y.^2)+1;
num=0;
num=num+1;
subplot(2,3,num);
plot3(X,Y,Z)
axis([-2 2 -2 2 0 2]); %限定显示的范围
xlabel('x轴');
ylabel('y轴');
zlabel('z轴');
title('练习')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num=num+1;  
subplot(2,3,num);  
mesh(X,Y,Z);  
axis([-3 3 -3 3 0 2]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(2)');%标题  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
num=num+1;  
subplot(2,3,num);  
meshc(X,Y,Z);  
axis([-3 3 -3 3 0 2]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(3)');%标题  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
num=num+1;  
subplot(2,3,num);  
surf(X,Y,Z);  
axis([-3 3 -3 3 0 2]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(4)');%标题  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
num=num+1;  
subplot(2,3,num);  
meshz(X,Y,Z);  
axis([-3 3 -3 3 0 2]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(5)');%标题  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
num=num+1;  
subplot(2,3,num);  
surf(X,Y,Z);  
hold on;  
stem3(X,Y,Z,'r');%画竖线  
axis([-3 3 -3 3 0 2]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(6)');%标题  

%%%%%% 画矩阵的三维图！！！
clc
clear all
X=[0 1 2 3 4 5 6 7 8 9];
Y=[0 1 2 3 4 5 6 7 8 9]; 
for i=1:1:length(X)  
    for j=1:1:length(Y)  
        Z(i,j)=mod(i*j*rand(1),9);  
    end  
end  
num=0;  
num=num+1;  
subplot(2,3,num);  
plot3(X,Y,Z);  
axis([0 9 0 9 0 9]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(1)');%标题

num=num+1;  
subplot(2,3,num);  
mesh(X,Y,Z);  
axis([0 9 0 9 0 9]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(2)');%标题  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
num=num+1;  
subplot(2,3,num);  
meshc(X,Y,Z);  
axis([0 9 0 9 0 9]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(3)');%标题  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
num=num+1;  
subplot(2,3,num);  
surf(X,Y,Z);  
axis([0 9 0 9 0 9]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(4)');%标题  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
num=num+1;  
subplot(2,3,num);  
meshz(X,Y,Z);  
axis([0 9 0 9 0 9]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(5)');%标题  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
num=num+1;  
subplot(2,3,num);  
surf(X,Y,Z);  
hold on;  
stem3(X,Y,Z,'r');%画竖线  
axis([0 9 0 9 0 9]);%限定显示的范围  
xlabel('x轴');%x轴坐标  
ylabel('y轴');%y轴坐标  
zlabel('z轴');%z轴坐标  
title('http://blog.csdn.net/nuptboyzhb/   figure(6)');%标题  
