function uh = Maximum_Limiter(uh)
global Nx NumGLP dimPk m M phiGLL

uhGLL = zeros(Nx,NumGLP);
for i = 1:Nx
    for d = 1:dimPk
        uhGLL(i,:) = uhGLL(i,:) + uh(i,d)*phiGLL(:,d)';
    end
end

% calculate theta for each cell I_i
for i = 1:Nx
    theta = 1;
    for i1 = 1:NumGLP
        thetaq = min(abs([ (M - uh(i,1))/(uhGLL(i,i1) - uh(i,1)), (m - uh(i,1))/(uhGLL(i,i1) - uh(i,1)), 1 ]));
        if thetaq < theta
            theta = thetaq;
        end
    end
    uh(i,2:end) = uh(i,2:end)*theta;
end

end