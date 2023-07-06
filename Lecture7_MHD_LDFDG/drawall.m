% drawall.m
figure(1);
mesh(rc,zc,Q1);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q1)),max(max(Q1)) + 1e-6]);
title('rho')
colormap(cool);

figure(2);
mesh(rc,zc,Q2);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q2)),max(max(Q2)) + 1e-6]);
title('ur')
colormap(cool);

figure(3);
mesh(rc,zc,Q3);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q3)),max(max(Q3)) + 1e-6]);
title('uz')
colormap(cool);

figure(4);
mesh(rc,zc,Q4);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q4)),max(max(Q4)) + 1e-6]);
title('uphi')
colormap(cool);

figure(5);
mesh(rc,zc,Q5);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q5)),max(max(Q5)) + 1e-6]);
title('E')
colormap(cool);

figure(6);
mesh(rc,zc,Q6);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q6)),max(max(Q6)) + 1e-6]);
title('Br')
colormap(cool);

figure(7);
mesh(rc,zc,Q7);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q7)),max(max(Q7)) + 1e-6]);
title('Bz')
colormap(cool);

figure(8);
mesh(rc,zc,Q8);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q8)),max(max(Q8)) + 1e-6]);
title('Bphi')
colormap(cool);