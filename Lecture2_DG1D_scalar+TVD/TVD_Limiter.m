function uh = TVD_Limiter(uh)

global Nx dimPk bcL bcR

uhb = [[0,0,0];uh;[0,0,0]];
uhmod = zeros(Nx,dimPk);
uhmod(:,1) = uh(:,1);

% set_bc
if bcL == 1
    uhb(1,:) = uh(end,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
end

for i = 1:Nx
    deltaUR = uh(i,2) + (2/3)*uh(i,3);
    deltaUL = uh(i,2) - (2/3)*uh(i,3);
    deltaURM = uhb(i + 2,1) - uhb(i + 1,1);
    deltaULM = uhb(i + 1,1) - uhb(i,1);
    
    deltaURM1 = minmod(deltaUR,deltaURM,deltaULM);
    deltaULM1 = minmod(deltaUL,deltaURM,deltaULM);
    
    uhmod(i,2) = (deltaURM1 + deltaULM1)/2;
    uhmod(i,3) = 3*(deltaURM1 - deltaULM1)/4;
end

uh = uhmod;

end

function a1 = minmod(a,b,c)

global hx M

if abs(a) < M*hx^2
    a1 = a;
else
    if sign(a) == sign(b) && sign(a) == sign(c)
        a1 = sign(a)*min(abs([a,b,c]));
    else
        a1 = 0;
    end
end

end