%average_to_GLpoint.m
if length(Xc) > 1
    hx1 = (Xc(2) - Xc(1))/2;
else
    hx1 = 1;
end
if length(Yc) > 1
    hy1 = (Yc(2) - Yc(1))/2;
else
    hy1 = 1;
end
NewGLP = 11;

NumGLP = NewGLP;

lambda = -1:0.2:1;

get_basis

XGflash = zeros(1,Nx*NumGLP);
YGflash = zeros(1,Ny*NumGLP);

for i = 1:Nx
    XGflash(NumGLP*(i - 1) + 1:NumGLP*i) = Xc(i) + hx1*lambda;
end

for j = 1:Ny
    YGflash(NumGLP*(j - 1) + 1:NumGLP*j) = Yc(j) + hy1*lambda;
end

xcG = zeros(Nx*NumGLP,Ny*NumGLP);
ycG = zeros(Nx*NumGLP,Ny*NumGLP);

for j = 1:Ny*NumGLP
    xcG(:,j) = XGflash;
end

for i = 1:Nx*NumGLP
    ycG(i,:) = YGflash;
end

Q1G = zeros(Nx*NumGLP,Ny*NumGLP);
Q2G = zeros(Nx*NumGLP,Ny*NumGLP);
Q3G = zeros(Nx*NumGLP,Ny*NumGLP);
Q4G = zeros(Nx*NumGLP,Ny*NumGLP);

for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:NumGLP
            for j1 = 1:NumGLP
                for d = 1:dimPk
                    Q1G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q1G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q1h(i,j,d)*phiG(i1,j1,d);
                    Q2G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q2G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q2h(i,j,d)*phiG(i1,j1,d);
                    Q3G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q3G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q3h(i,j,d)*phiG(i1,j1,d);
                    Q4G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q4G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q4h(i,j,d)*phiG(i1,j1,d);
                end
            end
        end
    end
end

Q1 = Q1G;
Q2 = Q2G;
Q3 = Q3G;
Q4 = Q4G;

Nx = Nx*NumGLP;
Ny = Ny*NumGLP;
Xc = XGflash;
Yc = YGflash;
xc = xcG;
yc = ycG;