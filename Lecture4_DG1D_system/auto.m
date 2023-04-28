% auto.m
clear;clc
global flux_type Nx Limit_type
Limit_type = 1;
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
Xc1 = Xc;

% calculate the "real" solution
Limit_type = 1;
Nx = 4000;
main
draw_solution1