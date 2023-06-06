% drawall.m
figure(1);
mesh(xc,yc,Q1);colorbar
title('rho')
colormap(cool);

figure(2);
mesh(xc,yc,Q2./Q1);
title('rhou')
colormap(cool);

figure(3);
mesh(xc,yc,Q3);
title('rhov')
colormap(cool);

figure(4);
mesh(xc,yc,Q4);
title('E')
colormap(cool);