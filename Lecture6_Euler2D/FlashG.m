% FlashG.m

h = figure();
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(h,'JavaFrame');
pause(0.1);
set(jFrame,'Maximized',1);
pause(0.1);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

QF = Q1flashG;
FRAME = frame;
t0 = T(end)/FRAME;

for i = 1:FRAME
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    mesh(xcG,ycG,QF(:,:,j));colormap(cool);axis([XcG(1),XcG(end),YcG(1),YcG(end),min(min(min(QF))) - 0.1,max(max(max(QF))) + 0.1]);
    %contour(xcG,ycG,QF(:,:,j),15);colormap(nclCM(484,100));
    %p = pcolor(xcG,ycG,QF(:,:,j));p.EdgeColor = 'none';colormap(nclCM(484,100));
    %imshow(QF(:,:,j)'/max(max(max(abs(QF(:,:,j))))))
    %colormap(jet)
    pause(0.0001);
end