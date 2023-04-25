function du = Lh1(uh)
global Nx bcL bcR hx NumEq flux_type
uhb = zeros(Nx + 2,NumEq);
uhb(2:end - 1,:,:) = uh;
uhR = zeros(Nx + 1,NumEq);
uhL = zeros(Nx + 1,NumEq);
fhat = zeros(Nx + 1,NumEq);
du = zeros(Nx,NumEq);

% set_bc
if bcL == 1
    uhb(1,:) = uh(end,:);
elseif bcL == 2
    uhb(1,:) = uh(1,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
elseif bcR == 2
    uhb(end,:) = uh(end,:);
end

% Step 2: calculate the flux at edge
for i = 1:Nx + 1
    for n = 1:NumEq
        uhR(i,n) = uhb(i,n);
        uhL(i,n) = uhb(i + 1,n);
    end
end

for i = 1:Nx + 1
    uR = uhL(i,:);
    uL = uhR(i,:);
    fR(1) = f1(uR(1),uR(2),uR(3));
    fR(2) = f2(uR(1),uR(2),uR(3));
    fR(3) = f3(uR(1),uR(2),uR(3));
    fL(1) = f1(uL(1),uL(2),uL(3));
    fL(2) = f2(uL(1),uL(2),uL(3));
    fL(3) = f3(uL(1),uL(2),uL(3));
    [SLmax,SLmin] = wavespeed(uL);
    [SRmax,SRmin] = wavespeed(uR);
    SR = max(SLmax,SRmax);
    SL = min(SLmin,SRmin);
    switch flux_type
        case 1
            fhat(i,:) = LF_Flux(uR,uL,fR,fL,SR,SL);
        case 2
            fhat(i,:) = HLL_Flux(uR,uL,fR,fL,SR,SL);
        case 3
            fhat(i,:) = HLLC_Flux(uR,uL,fR,fL,SR,SL);
    end
end

% Step 3: calculate [f*phi] on each cell
for i = 1:Nx
    for n = 1:NumEq
        du(i,n) = du(i,n) - (1/hx)*(fhat(i + 1,n) - fhat(i,n));
    end
end

end