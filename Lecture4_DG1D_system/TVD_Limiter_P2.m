function uh = TVD_Limiter_P2(uh)

global Nx dimPk bcL bcR NumEq gamma Limit_type

uhb = zeros(Nx + 2,dimPk,NumEq);
uhb(2:end - 1,:,:) = uh;
uhmod = zeros(Nx,dimPk,NumEq);
uhmod(:,1,:) = uh(:,1,:);
deltaUR = zeros(NumEq,1);
deltaUL = zeros(NumEq,1);
deltaURM = zeros(NumEq,1);
deltaULM = zeros(NumEq,1);

% set_bc
if bcL == 1
    uhb(1,:,:) = uh(end,:,:);
elseif bcL == 2
    uhb(1,:,:) = uh(2,:,:);
end

if bcR == 1
    uhb(end,:,:) = uh(1,:,:);
elseif bcR == 2
    uhb(end,:,:) = uh(end - 1,:,:);
end

for i = 1:Nx
    for n = 1:NumEq
        deltaUR(n,1) = uh(i,2,n) + (2/3)*uh(i,3,n);
        deltaUL(n,1) = uh(i,2,n) - (2/3)*uh(i,3,n);
        deltaURM(n,1) = uhb(i + 2,1,n) - uhb(i + 1,1,n);
        deltaULM(n,1) = uhb(i + 1,1,n) - uhb(i,1,n);
    end
    
    v = uh(i,1,2)/uh(i,1,1);
    p = pressure(uh(i,1,1),uh(i,1,2),uh(i,1,3));
    c = sqrt(gamma*p/uh(i,1,1));
    H = (uh(i,1,3) + p)/uh(i,1,1);
    
    R1 = [1;v - c;H - v*c];
    R2 = [1;v;0.5*v^2];
    R3 = [1;v + c;H + v*c];
    
    R = [R1,R2,R3];%+ 1e-12*eye(3);
     
    L = inv(R);
    
    deltaURM1 = R*minmod(L*deltaUR,L*deltaURM,L*deltaULM);
    deltaULM1 = R*minmod(L*deltaUL,L*deltaURM,L*deltaULM);
    
    for n = 1:NumEq
        uhmod(i,2,n) = (deltaURM1(n) + deltaULM1(n))/2;
        uhmod(i,3,n) = 3*(deltaURM1(n) - deltaULM1(n))/4;
    end
end

uh = uhmod;

end

function a1 = minmod(a,b,c)

global hx M NumEq

a1 = zeros(NumEq,1);

for i = 1:NumEq
    if abs(a(i)) < M*hx^2
        a1(i) = a(i);
    else
        if sign(a(i)) == sign(b(i)) && sign(a(i)) == sign(c(i))
            a1(i) = sign(a(i))*min(abs([a(i),b(i),c(i)]));
        else
            a1(i) = 0;
        end
    end
end

end