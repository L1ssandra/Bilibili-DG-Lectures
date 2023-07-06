% drawMachstem.m
% x = 0.3 to 0.7, and y = 0.3 to 0.7
Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
QE = Q4;
QB1 = Q5;
QB2 = Q6;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP = 0.5*(QB1.^2 + QB2.^2);

QMach = QMach(:,:,end);
Nx = length(Xc); Ny = length(Yc);
XcM = [];
YcM = [];
i1 = Nx; i2 = 1;
j1 = Ny; j2 = 1;
for i = 1:Nx
    if Xc(i) > 0.3 && Xc(i) < 0.7
        XcM = [XcM,X(i)];
        i1 = min([i,i1]);
        i2 = max([i,i2]);
    end
end

for j = 1:Ny
    if Yc(j) > 0.3 && Yc(j) < 0.7
        YcM = [YcM,Y(j)];
        j1 = min([j,j1]);
        j2 = max([j,j2]);
    end
end

QMachM = QMach(i1:i2,j1:j2);
NxM = length(XcM);
NyM = length(YcM);

xcM = zeros(NxM,NyM);
ycM = zeros(NxM,NyM);

for j = 1:NyM
    xcM(:,j) = XcM;
end

for i = 1:NxM
    ycM(i,:) = YcM';
end

figure(1);
contour(xcM,ycM,QMachM,30);colormap(cool);colorbar
%p = pcolor(xcM,ycM,Q1M);p.EdgeColor = 'none';colormap(jet);colorbar;
title('Mach number')