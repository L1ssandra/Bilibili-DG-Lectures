% drawlogcontour.m
figure(1);
p = pcolor(xc,yc,log(Q1));p.EdgeColor = 'none';colormap(jet);colorbar;
%contour(xc,yc,Q1,40);colormap(cool)
title('density')

figure(2);
gamma = 1.4;
p = pcolor(xc,yc,log((gamma - 1)*(Q4 - 0.5*(Q2.^2 + Q3.^2)./Q1)));p.EdgeColor = 'none';colormap(jet);colorbar;
%contour(xc,yc,Q1,40);colormap(cool)
title('pressure')