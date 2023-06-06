% drawMachstem.m
% x = 2.2 to 2.8, and y = 0 to 0.5
XcM = [];
YcM = [];
i1 = Nx; i2 = 1;
j1 = 1; j2 = 1;
for i = 1:Nx
    if Xc(i) > 2.2 && Xc(i) < 2.8
        XcM = [XcM,Xc(i)];
        i1 = min([i,i1]);
        i2 = max([i,i2]);
    end
end

for j = 1:Ny
    if Yc(j) < 0.5
        YcM = [YcM,Yc(j)];
        j2 = max([j,j2]);
    end
end

Q1M = Q1(i1:i2,j1:j2);
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
%contourf(xcM,ycM,Q1M,20);colormap(nclCM(203,150));colorbar
p = pcolor(xcM,ycM,Q1M);p.EdgeColor = 'none';colormap(nclCM(203,150));colorbar;
title('density')