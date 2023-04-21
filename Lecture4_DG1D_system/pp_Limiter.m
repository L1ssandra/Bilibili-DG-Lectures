function uh = pp_Limiter(uh)
global Nx NumGLP dimPk phiGLL NumEq

epsilon = 1e-12;

uhGLL = zeros(Nx,NumGLP,NumEq);
for i = 1:Nx
    for d = 1:dimPk
        for n = 1:NumEq
            for i1 = 1:NumGLP
                uhGLL(i,i1,n) = uhGLL(i,i1,n) + uh(i,d,n)*phiGLL(i1,d);
            end
        end
    end
end

% Limiting the density
for i = 1:Nx
    theta = 1;
    for i1 = 1:NumGLP
        thetaq = min(abs([ (uh(i,1,1) - epsilon)/(uh(i,1,1) - uhGLL(i,i1,1)), 1 ]));
        if thetaq < theta
            theta = thetaq;
        end
    end
    uh(i,2:end,1) = uh(i,2:end,1)*theta;
end

uhGLL = zeros(Nx,NumGLP,NumEq);
for i = 1:Nx
    for d = 1:dimPk
        for n = 1:NumEq
            for i1 = 1:NumGLP
                uhGLL(i,i1,n) = uhGLL(i,i1,n) + uh(i,d,n)*phiGLL(i1,d);
            end
        end
    end
end

% Limiting the pressure
for i = 1:Nx
    theta = 1;
    for i1 = 1:NumGLP
        pq = pressure(uhGLL(i,i1,1),uhGLL(i,i1,2),uhGLL(i,i1,3));
        if pq < epsilon
            pt = @(t) pressure((1 - t)*uh(i,1,1) + t*uhGLL(i,i1,1),(1 - t)*uh(i,1,2) + t*uhGLL(i,i1,2),(1 - t)*uh(i,1,3) + t*uhGLL(i,i1,3)) - epsilon;
            tq = bisect(pt,0,1,epsilon);
            if tq < theta
                theta = tq;
            end
        end
    end
    
    if theta < 1
        theta = 0.95*theta;
    end
    uh(i,2:end,:) = uh(i,2:end,:)*theta;
end

end