% init_data.m
global bcL bcR hx hx1 M gamma

%----------------
xa = 0;
xb = 1;

gamma = 1.4;

rho = @(x) rho0(x);%1 + 0.2*sin(x);
u = @(x) u0(x);%1 + 0.*x;
p = @(x) p0(x);%1 + 0.*x;

U1 = @(x) rho(x); % rho
U2 = @(x) rho(x).*u(x); % rho u
U3 = @(x) p(x)./(gamma - 1) + 0.5*rho(x).*u(x).^2; % E

M = 1;
bcL = 2;
bcR = 2;
tend = 0.012;
%----------------

hx = (xb - xa)/Nx;
hx1 = 0.5*hx;

Xc = xa + hx1:hx:xb - hx1;

ureal = zeros(Nx,NumGLP,NumEq);
for i = 1:Nx
    for j = 1:NumGLP
        ureal(i,j,1) = U1(Xc(i) + hx1*lambda(j));
        ureal(i,j,2) = U2(Xc(i) + hx1*lambda(j));
        ureal(i,j,3) = U3(Xc(i) + hx1*lambda(j));
    end
end