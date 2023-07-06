% contourall.m
figure(1);
contour(rc,zc,Q1,20);colormap(cool);colorbar
title('rho')
colormap(cool);

figure(2);
contour(rc,zc,Q2,20);colormap(cool);colorbar
title('ur')
colormap(cool);

figure(3);
contour(rc,zc,Q3,20);colormap(cool);colorbar
title('uz')
colormap(cool);

figure(4);
contour(rc,zc,Q4,20);colormap(cool);colorbar
title('uphi')
colormap(cool);

figure(5);
contour(rc,zc,Q5,20);colormap(cool);colorbar
title('E')
colormap(cool);

figure(6);
contour(rc,zc,Q6,20);colormap(cool);colorbar
title('Br')
colormap(cool);

figure(7);
contour(rc,zc,Q7,20);colormap(cool);colorbar
title('Bz')
colormap(cool);

figure(8);
contour(rc,zc,Q8,20);colormap(cool);colorbar
title('Bphi')
colormap(cool);