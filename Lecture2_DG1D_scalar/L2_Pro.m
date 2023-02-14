% L2_Pro.m
uh = zeros(Nx,dimPk);

for i = 1:Nx
    for d = 1:dimPk
        for i1 = 1:NumGLP
            uh(i,d) = uh(i,d) + 0.5*weight(i1)*ureal(i,i1)*phiG(i1,d);
        end
    end
end

for d = 1:dimPk
    uh(:,d) = uh(:,d)/mm(d);
end