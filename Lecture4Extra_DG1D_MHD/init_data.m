% init_data.m
global bcL bcR hx hx1 M gamma

%----------------
xa = -1;
xb = 1;

gamma = 5/3;

rho = @(x) rho0(x);
u = @(x) 0.*x;
v = @(x) 0.*x;
p = @(x) p0(x);
Bx = @(x) 0.75 + 0.*x;
By = @(x) By0(x);

U1 = @(x) rho(x); % rho
U2 = @(x) rho(x).*u(x); % rho u
U3 = @(x) rho(x).*v(x); % rho v
U4 = @(x) p(x)./(gamma - 1) + 0.5*rho(x).*(u(x).^2 + v(x).^2) + 0.5*(Bx(x).^2 + By(x).^2); % E
U5 = @(x) Bx(x); % Bx
U6 = @(x) By(x); % By

M = 1;
bcL = 2;
bcR = 2;
tend = 0.2;
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
        ureal(i,j,4) = U4(Xc(i) + hx1*lambda(j));
        ureal(i,j,5) = U5(Xc(i) + hx1*lambda(j));
        ureal(i,j,6) = U6(Xc(i) + hx1*lambda(j));
    end
end