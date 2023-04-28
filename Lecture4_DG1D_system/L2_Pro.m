% L2_Pro.m
uh = zeros(Nx,dimPk,NumEq);

for i = 1:Nx
    for d = 1:dimPk
        for n = 1:NumEq
            for i1 = 1:NumGLP
                uh(i,d,n) = uh(i,d,n) + 0.5*weight(i1)*ureal(i,i1,n)*phiG(i1,d);
            end
        end
    end
end

for d = 1:dimPk
    uh(:,d,:) = uh(:,d,:)/mm(d);
end