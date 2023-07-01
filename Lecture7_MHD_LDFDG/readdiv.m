% readdiv.m
RcG = load('Xc.txt');
ZcG = load('Yc.txt');
lambda = load('lambda.txt');
weight = load('weight.txt');
NumGLP = length(lambda);
NrG = length(RcG);
NzG = length(ZcG);
Nr = NrG/NumGLP;
Nz = NzG/NumGLP;
dimPk = 6;
Rc = zeros(1,Nr);
Zc = zeros(1,Nz);

for i = 1:Nr
    Rc(i) = (RcG((i - 1)*NumGLP + 1) + RcG(i*NumGLP))/2;
end

for j = 1:Nz
    Zc(j) = (ZcG((j - 1)*NumGLP + 1) + ZcG(j*NumGLP))/2;
end

Qdiv = load('div.txt');

rc = zeros(Nr,Nz);
zc = zeros(Nr,Nz);

for j = 1:Nz
    rc(:,j) = Rc;
end

for i = 1:Nr
    zc(i,:) = Zc';
end

Qdiv = reshape(Qdiv,Nr,Nz);

%mesh(rc,zc,Qdiv);axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Qdiv)),max(max(Qdiv))]);
p = pcolor(rc,zc,Qdiv); colormap(jet); p.EdgeColor = 'none'; colorbar
%imshow(Qdiv/max(max(Qdiv)))