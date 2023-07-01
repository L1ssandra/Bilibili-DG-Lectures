% Flash.m

h = figure(1);
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(h,'JavaFrame');
pause(0.1);
set(jFrame,'Maximized',1);
pause(0.1);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

Qrhoflash = Q1flash;
Quflash = Q2flash./Q1flash;
Qvflash = Q3flash./Q1flash;
Qwflash = Q4flash./Q1flash;
QEflash = Q5flash;
QB1flash = Q6flash;
QB2flash = Q7flash;
QB3flash = Q8flash;
gamma = 5/3;
QPflash = (gamma - 1)*(QEflash - 0.5*Qrhoflash.*(Quflash.^2 + Qvflash.^2 + Qwflash.^2) - 0.5*(QB1flash.^2 + QB2flash.^2 + QB3flash.^2));
QCflash = sqrt(abs(gamma*QPflash./Qrhoflash));
QMachflash = sqrt(Quflash.^2 + Qvflash.^2 + Qwflash.^2)./QCflash;
QBPflash = 0.5*(QB1flash.^2 + QB2flash.^2 + QB3flash.^2);
QV2flash = Quflash.^2 + Qvflash.^2 + Qwflash.^2;

QF = QBPflash;
FRAME = frame;
t0 = T(end)/(FRAME - 1);

for i = 1:FRAME
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    %contourf(rc,zc,QF(:,:,1,j),15);colormap(nclCM(484,100));
    %mesh(rc,zc,QF(:,:,1,j));axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(min(min(QF)))) - 0.1,max(max(max(max(QF)))) + 0.1]);colormap(cool);
    %mesh(rc,zc,QF(:,:,1,j));axis([Rc(1),Rc(end),Zc(1),Zc(end),0,2]);colormap(cool);
    p = pcolor(rc,zc,QF(:,:,1,j));colormap(nclCM(484,100));p.EdgeColor = 'none';
    colorbar;
    pause(0.0001);
end