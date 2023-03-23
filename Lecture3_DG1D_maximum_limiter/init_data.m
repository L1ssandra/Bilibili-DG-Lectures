% init_data.m
global bcL bcR hx hx1 m M

%----------------
xa = -1;
xb = 1;
%u0 = @(x) sin(x);
m = 0; M = 1; Mtvb = 0;
bcL = 1;
bcR = 1;
tend = 0.4;
%----------------

hx = (xb - xa)/Nx;
hx1 = 0.5*hx;

Xc = xa + hx1:hx:xb - hx1;

ureal = zeros(Nx,NumGLP);
for i = 1:Nx
    for j = 1:NumGLP
        ureal(i,j) = u0(Xc(i) + hx1*lambda(j));
    end
end