%FlashFortran.m
RcG = load('Rc.txt');
ZcG = load('Zc.txt');
lambda = load('lambda.txt');
weight = load('weight.txt');
NumGLP = length(lambda);
NrG = length(RcG);
NzG = length(ZcG);
Nr = NrG/NumGLP;
Nz = NzG/NumGLP;
dimPk = 10;
Rc = zeros(1,Nr);
Zc = zeros(1,Nz);

for i = 1:Nr
    Rc(i) = (RcG((i - 1)*NumGLP + 1) + RcG(i*NumGLP))/2;
end

for j = 1:Nz
    Zc(j) = (ZcG((j - 1)*NumGLP + 1) + ZcG(j*NumGLP))/2;
end

Q1flash = load('Q1flash.txt');
Q2flash = load('Q2flash.txt');
Q3flash = load('Q3flash.txt');
Q4flash = load('Q4flash.txt');
Q5flash = load('Q5flash.txt');
Q6flash = load('Q6flash.txt');
Q7flash = load('Q7flash.txt');
Q8flash = load('Q8flash.txt');
T = load('T.txt');

frame = length(T);

rc = zeros(Nr,Nz);
zc = zeros(Nr,Nz);

for j = 1:Nz
    rc(:,j) = Rc;
end

for i = 1:Nr
    zc(i,:) = Zc';
end

Q1hflash = reshape(Q1flash,Nr,Nz,dimPk,frame);
Q2hflash = reshape(Q2flash,Nr,Nz,dimPk,frame);
Q3hflash = reshape(Q3flash,Nr,Nz,dimPk,frame);
Q4hflash = reshape(Q4flash,Nr,Nz,dimPk,frame);
Q5hflash = reshape(Q5flash,Nr,Nz,dimPk,frame);
Q6hflash = reshape(Q6flash,Nr,Nz,dimPk,frame);
Q7hflash = reshape(Q7flash,Nr,Nz,dimPk,frame);
Q8hflash = reshape(Q8flash,Nr,Nz,dimPk,frame);

Q1flash = Q1hflash(:,:,1,:);
Q2flash = Q2hflash(:,:,1,:);
Q3flash = Q3hflash(:,:,1,:);
Q4flash = Q4hflash(:,:,1,:);
Q5flash = Q5hflash(:,:,1,:);
Q6flash = Q6hflash(:,:,1,:);
Q7flash = Q7hflash(:,:,1,:);
Q8flash = Q8hflash(:,:,1,:);