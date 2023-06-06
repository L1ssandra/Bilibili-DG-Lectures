% get_basis.m



phiG = zeros(NumGLP,NumGLP,dimPk);
phiGR = zeros(NumGLP,dimPk);
phiGL = zeros(NumGLP,dimPk);
phiGU = zeros(NumGLP,dimPk);
phiGD = zeros(NumGLP,dimPk);

for i = 1:NumGLP
    for j = 1:NumGLP
        phiG(i,j,1) = 1;
        phiGR(j,1) = 1;
        phiGL(j,1) = 1;
        phiGU(i,1) = 1;
        phiGD(i,1) = 1;
        
        phiG(i,j,2) = lambda(i);
        phiGR(j,2) = 1;
        phiGL(j,2) = -1;
        phiGU(i,2) = lambda(i);
        phiGD(i,2) = lambda(i);
        
        phiG(i,j,3) = lambda(j);
        phiGR(j,3) = lambda(j);
        phiGL(j,3) = lambda(j);
        phiGU(i,3) = 1;
        phiGD(i,3) = -1;
        
        phiG(i,j,4) = lambda(i)^2 - 1/3;
        phiGR(j,4) = 2/3;
        phiGL(j,4) = 2/3;
        phiGU(i,4) = lambda(i)^2 - 1/3;
        phiGD(i,4) = lambda(i)^2 - 1/3;
        
        phiG(i,j,5) = lambda(i)*lambda(j);
        phiGR(j,5) = lambda(j);
        phiGL(j,5) = -lambda(j);
        phiGU(i,5) = lambda(i);
        phiGD(i,5) = -lambda(i);
        
        phiG(i,j,6) = lambda(j)^2 - 1/3;
        phiGR(j,6) = lambda(j)^2 - 1/3;
        phiGL(j,6) = lambda(j)^2 - 1/3;
        phiGU(i,6) = 2/3;
        phiGD(i,6) = 2/3;
        
        phiG(i,j,7) = lambda(i)^3 - 3*lambda(i)/5;
        phiGR(j,7) = 2/5;
        phiGL(j,7) = -2/5;
        phiGU(i,7) = lambda(i)^3 - 3*lambda(i)/5;
        phiGD(i,7) = lambda(i)^3 - 3*lambda(i)/5;
        
        phiG(i,j,8) = (lambda(i)^2 - 1/3)*(lambda(j));
        phiGR(j,8) = (2/3)*(lambda(j));
        phiGL(j,8) = (2/3)*(lambda(j));
        phiGU(i,8) = (lambda(i)^2 - 1/3);
        phiGD(i,8) = -(lambda(i)^2 - 1/3);
        
        phiG(i,j,9) = (lambda(i))*(lambda(j)^2 - 1/3);
        phiGR(j,9) = (lambda(j)^2 - 1/3);
        phiGL(j,9) = -(lambda(j)^2 - 1/3);
        phiGU(i,9) = lambda(i)*(2/3);
        phiGD(i,9) = lambda(i)*(2/3);
        
        phiG(i,j,10) = lambda(j)^3 - 3*lambda(j)/5;
        phiGR(j,10) = lambda(j)^3 - 3*lambda(j)/5;
        phiGL(j,10) = lambda(j)^3 - 3*lambda(j)/5;
        phiGU(i,10) = 2/5;
        phiGD(i,10) = -2/5;
    end
end