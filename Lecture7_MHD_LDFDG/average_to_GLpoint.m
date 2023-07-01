%average_to_GLpoint.m
hr1 = (Rc(2) - Rc(1))/2;
hz1 = (Zc(2) - Zc(1))/2;
NewGLP = 11;

NumGLP = NewGLP;

lambda = -1:0.2:1;

get_basis

XGflash = zeros(1,Nr*NumGLP);
YGflash = zeros(1,Nz*NumGLP);

for i = 1:Nr
    XGflash(NumGLP*(i - 1) + 1:NumGLP*i) = Rc(i) + hr1*lambda;
end

for j = 1:Nz
    YGflash(NumGLP*(j - 1) + 1:NumGLP*j) = Zc(j) + hz1*lambda;
end

xcG = zeros(Nr*NumGLP,Nz*NumGLP);
ycG = zeros(Nr*NumGLP,Nz*NumGLP);

for j = 1:Nz*NumGLP
    xcG(:,j) = XGflash;
end

for i = 1:Nr*NumGLP
    ycG(i,:) = YGflash;
end

Q1G = zeros(Nr*NumGLP,Nz*NumGLP);
Q2G = zeros(Nr*NumGLP,Nz*NumGLP);
Q3G = zeros(Nr*NumGLP,Nz*NumGLP);
Q4G = zeros(Nr*NumGLP,Nz*NumGLP);
Q5G = zeros(Nr*NumGLP,Nz*NumGLP);
Q6G = zeros(Nr*NumGLP,Nz*NumGLP);
Q7G = zeros(Nr*NumGLP,Nz*NumGLP);
Q8G = zeros(Nr*NumGLP,Nz*NumGLP);

for i = 1:Nr
    for j = 1:Nz
        for i1 = 1:NumGLP
            for j1 = 1:NumGLP
                for d = 1:dimPk
                    Q1G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q1G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q1h(i,j,d)*phiG(i1,j1,d);
                    Q2G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q2G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q2h(i,j,d)*phiG(i1,j1,d);
                    Q3G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q3G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q3h(i,j,d)*phiG(i1,j1,d);
                    Q4G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q4G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q4h(i,j,d)*phiG(i1,j1,d);
                    Q5G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q5G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q5h(i,j,d)*phiG(i1,j1,d);
                    Q6G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q6G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q6h(i,j,d)*phiG(i1,j1,d);
                    Q7G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q7G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q7h(i,j,d)*phiG(i1,j1,d);
                    Q8G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q8G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q8h(i,j,d)*phiG(i1,j1,d);
                end
            end
        end
    end
end

Q1 = Q1G;
Q2 = Q2G;
Q3 = Q3G;
Q4 = Q4G;
Q5 = Q5G;
Q6 = Q6G;
Q7 = Q7G;
Q8 = Q8G;

Rc = XGflash;
Zc = YGflash;
rc = xcG;
zc = ycG;