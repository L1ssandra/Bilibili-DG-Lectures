% drawcontourf.m
figure(1);
contourf(xc,yc,Q1,20);colormap(nclCM(203,150));colorbar
title('density')

% figure(2);
% gamma = 1.4;
% contourf(xc,yc,Q1,12);colormap(jet);colorbar
% title('pressure')