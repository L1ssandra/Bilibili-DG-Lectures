function du = Lh(uh)
global Nx NumGLP dimPk mm phiG phixG phiGR phiGL bcL bcR weight hx
uhb = [[0,0,0];uh;[0,0,0]];
uhG = zeros(Nx,NumGLP);
uhR = zeros(Nx + 1,1);
uhL = zeros(Nx + 1,1);
fhat = zeros(Nx + 1,1);
du = zeros(Nx,dimPk);

% set_bc
if bcL == 1
    uhb(1,:) = uh(end,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
end

% Step 1: calculate the Integral in cell
for i = 1:Nx
    for d = 1:dimPk
        uhG(i,:) = uhG(i,:) + uh(i,d)*phiG(:,d)';
    end
end

for i = 1:Nx
    for d = 2:dimPk
        for i1 = 1:NumGLP
            du(i,d) = du(i,d) + 0.5*weight(i1)*f(uhG(i,i1))*phixG(i1,d);
        end
    end
end

% Step 2: calculate the flux at edge
for i = 1:Nx + 1
    for d = 1:dimPk
        uhR(i) = uhR(i) + uhb(i,d)*phiGR(d);
        uhL(i) = uhL(i) + uhb(i + 1,d)*phiGL(d);
    end
end

for i = 1:Nx + 1
    uR = uhL(i);
    uL = uhR(i);
    alpha = 1;
    fhat(i) = 0.5*(f(uR) + f(uL) - alpha*(uR - uL));
end

% Step 3: calculate [f*phi] on each cell
for i = 1:Nx
    for d = 1:dimPk
        du(i,d) = du(i,d) - (1/hx)*(phiGR(d)*fhat(i + 1) - phiGL(d)*fhat(i));
    end
end

for d = 1:dimPk
    du(:,d) = du(:,d)/mm(d);
end

