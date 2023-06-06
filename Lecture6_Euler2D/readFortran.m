% readFortran.m

XcG = load('Xc.txt');
YcG = load('Yc.txt');
lambda = load('lambda.txt');
weight = load('weight.txt');
NumGLP = length(lambda);
NxG = length(XcG);
NyG = length(YcG);
Nx = NxG/NumGLP;
Ny = NyG/NumGLP;
dimPk = 6;
Xc = zeros(1,Nx);
Yc = zeros(1,Ny);

for i = 1:Nx
    Xc(i) = (XcG((i - 1)*NumGLP + 1) + XcG(i*NumGLP))/2;
end

for j = 1:Ny
    Yc(j) = (YcG((j - 1)*NumGLP + 1) + YcG(j*NumGLP))/2;
end

Q1 = load('Q1.txt');
Q2 = load('Q2.txt');
Q3 = load('Q3.txt');
Q4 = load('Q4.txt');

xc = zeros(Nx,Ny);
yc = zeros(Nx,Ny);

for j = 1:Ny
    xc(:,j) = Xc;
end

for i = 1:Nx
    yc(i,:) = Yc';
end

Q1h = reshape(Q1,Nx,Ny,dimPk);
Q2h = reshape(Q2,Nx,Ny,dimPk);
Q3h = reshape(Q3,Nx,Ny,dimPk);
Q4h = reshape(Q4,Nx,Ny,dimPk);

%drawRotor
%drawSmoothVortex
%drawOTV
%drawall
%drawdiv
Q1 = Q1h(:,:,1);
Q2 = Q2h(:,:,1);
Q3 = Q3h(:,:,1);
Q4 = Q4h(:,:,1);
%draw1