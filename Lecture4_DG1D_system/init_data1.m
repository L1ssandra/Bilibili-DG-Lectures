global bcL bcR hx hx1 gamma

%----------------
xa = 0;
xb = 1;

gamma = 1.4;

rho = @(x) rho0(x);
u = @(x) u0(x);
p = @(x) p0(x);

U1 = @(x) rho(x); % rho
U2 = @(x) rho(x).*u(x); % rho u
U3 = @(x) p(x)./(gamma - 1) + 0.5*rho(x).*u(x).^2; % E

bcL = 2;
bcR = 2;
tend = 0.012;
%----------------

hx = (xb - xa)/Nx;
hx1 = 0.5*hx;

Xc = xa + hx1:hx:xb - hx1;

uh = zeros(Nx,NumEq);

for i = 1:Nx
    for n = 1:NumEq
        uh(i,1) = U1(Xc(i));
        uh(i,2) = U2(Xc(i));
        uh(i,3) = U3(Xc(i));
    end
end