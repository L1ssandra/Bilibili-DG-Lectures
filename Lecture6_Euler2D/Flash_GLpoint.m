% Flash_GLpoint
hx1 = (Xc(2) - Xc(1))/2;
hy1 = (Yc(2) - Yc(1))/2;
NewGLP = 5;

NumGLP = NewGLP;

if NewGLP == 10
    lambda = [0.1488743389816312108848260;
        0.4333953941292471907992659;
        0.6794095682990244062343274;
        0.8650633666889845107320967;
        0.9739065285171717200779640;
        -0.1488743389816312108848260;
        -0.4333953941292471907992659;
        -0.6794095682990244062343274;
        -0.8650633666889845107320967;
        -0.9739065285171717200779640];
elseif NewGLP == 5
    lambda = [-0.9061798459386639927976269;
        -0.5384693101056830910363144;
        0;
        0.5384693101056830910363144;
        0.9061798459386639927976269];
end

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

Q1flashG = zeros(Nx*NumGLP,Ny*NumGLP,frame);
Q2flashG = zeros(Nx*NumGLP,Ny*NumGLP,frame);
Q3flashG = zeros(Nx*NumGLP,Ny*NumGLP,frame);
Q4flashG = zeros(Nx*NumGLP,Ny*NumGLP,frame);
Q5flashG = zeros(Nx*NumGLP,Ny*NumGLP,frame);
Q6flashG = zeros(Nx*NumGLP,Ny*NumGLP,frame);
Q7flashG = zeros(Nx*NumGLP,Ny*NumGLP,frame);
Q8flashG = zeros(Nx*NumGLP,Ny*NumGLP,frame);

for t = 1:frame
    for i = 1:Nx
        for j = 1:Ny
            for i1 = 1:NumGLP
                for j1 = 1:NumGLP
                    for d = 1:dimPk
                        Q1flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) = Q1flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) + Q1hflash(i,j,d,t)*phiG(i1,j1,d);
                        Q2flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) = Q2flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) + Q2hflash(i,j,d,t)*phiG(i1,j1,d);
                        Q3flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) = Q3flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) + Q3hflash(i,j,d,t)*phiG(i1,j1,d);
                        Q4flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) = Q4flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) + Q4hflash(i,j,d,t)*phiG(i1,j1,d);
                        Q5flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) = Q5flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) + Q5hflash(i,j,d,t)*phiG(i1,j1,d);
                        Q6flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) = Q6flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) + Q6hflash(i,j,d,t)*phiG(i1,j1,d);
                        Q7flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) = Q7flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) + Q7hflash(i,j,d,t)*phiG(i1,j1,d);
                        Q8flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) = Q8flashG((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1,t) + Q8hflash(i,j,d,t)*phiG(i1,j1,d);
                    end
                end
            end
        end
    end
end

QrhoflashG = Q1flashG;
QuflashG = Q2flashG./Q1flashG;
QvflashG = Q3flashG./Q1flashG;
QwflashG = Q4flashG./Q1flashG;
QEflashG = Q5flashG;
QB1flashG = Q6flashG;
QB2flashG = Q7flashG;
QB3flashG = Q8flashG;
gamma = 5/3;
QPflashG = (gamma - 1)*(QEflashG - 0.5*QrhoflashG.*(QuflashG.^2 + QvflashG.^2 + QwflashG.^2) - 0.5*(QB1flashG.^2 + QB2flashG.^2 + QB3flashG.^2));
QCflashG = sqrt(abs(gamma*QPflashG./QrhoflashG));
QMachflashG = sqrt(QuflashG.^2 + QvflashG.^2 + QwflashG.^2)./QCflashG;
QBPflashG = 0.5*(QB1flashG.^2 + QB2flashG.^2 + QB3flashG.^2);
QBnormflashG = (2*QBPflashG).^0.5;
QV2flashG = QuflashG.^2 + QvflashG.^2 + QwflashG.^2;