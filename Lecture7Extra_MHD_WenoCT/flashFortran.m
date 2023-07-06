% 显示解随时间变化规律,从Fortran写好的文件中读取数据
% flashFortran.m
clear;clc
xa = 0; xb = 1; ya = 0; yb = 1;

Q1 = load('Q1.txt');
Q2 = load('Q2.txt');
Q3 = load('Q3.txt');
Q4 = load('Q4.txt');
Q5 = load('Q5.txt');
Q6 = load('Q6.txt');

Q1 = Q1(:,1);
Q2 = Q2(:,1);
Q3 = Q3(:,1);
Q4 = Q4(:,1);
Q5 = Q5(:,1);
Q6 = Q6(:,1);

T = load('T.txt');

Nx = load('Nx.txt');
Ny = load('Ny.txt');

frame = length(T);

N = round(sqrt(length(Q1)/frame));
% 网格
hx = (xb - xa)/Nx;
hy = (yb - ya)/Ny;

X = zeros(Nx,1);
Y = zeros(Ny,1);
for i = 1:Nx + 1
    X(i) = xa + (i - 1)*hx;
end

for i = 1:Ny + 1
    Y(i) = ya + (i - 1)*hy;
end

% 整格点
Xc = (X(1:end - 1) + X(2:end))/2;
Yc = (Y(1:end - 1) + Y(2:end))/2;

Q1 = reshape(Q1,Nx,Ny,frame);
Q2 = reshape(Q2,Nx,Ny,frame);
Q3 = reshape(Q3,Nx,Ny,frame);
Q4 = reshape(Q4,Nx,Ny,frame);
Q5 = reshape(Q5,Nx,Ny,frame);
Q6 = reshape(Q6,Nx,Ny,frame);

[xc,yc] = meshgrid(Yc,Xc);
%plot(T,DIV);
% contour(Xc,Yc,QF(:,:,j,1),25);
%draw4
%flash2