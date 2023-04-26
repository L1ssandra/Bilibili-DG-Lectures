% auto.m
global flux_type Nx
Nx = 200;
flux_type = 1;
main
u1 = uh;
flux_type = 2;
main
u2 = uh;
flux_type = 3;
main
u3 = uh;
flux_type = 4;
main
u4 = uh;
Xc1 = Xc;

Nx = 1000;
main
draw_solution1
