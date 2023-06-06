%FlashFortran.m
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

Q1flash = load('Q1flash.txt');
Q2flash = load('Q2flash.txt');
Q3flash = load('Q3flash.txt');
Q4flash = load('Q4flash.txt');
T = load('T.txt');

frame = length(T);

xc = zeros(Nx,Ny);
yc = zeros(Nx,Ny);

for j = 1:Ny
    xc(:,j) = Xc;
end

for i = 1:Nx
    yc(i,:) = Yc';
end

Q1hflash = reshape(Q1flash,Nx,Ny,dimPk,frame);
Q2hflash = reshape(Q2flash,Nx,Ny,dimPk,frame);
Q3hflash = reshape(Q3flash,Nx,Ny,dimPk,frame);
Q4hflash = reshape(Q4flash,Nx,Ny,dimPk,frame);

Q1flash = Q1hflash(:,:,1,:);
Q2flash = Q2hflash(:,:,1,:);
Q3flash = Q3hflash(:,:,1,:);
Q4flash = Q4hflash(:,:,1,:);

Qrhoflash = Q1flash;
Quflash = Q2flash./Q1flash;
Qvflash = Q3flash./Q1flash;
QEflash = Q4flash;
gamma = 5/3;
QPflash = (gamma - 1)*(QEflash - 0.5*Qrhoflash.*(Quflash.^2 + Qvflash.^2 ));
QCflash = sqrt(abs(gamma*QPflash./Qrhoflash));