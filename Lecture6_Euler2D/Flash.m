% Flash.m

% h = figure(1);
% warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
% jFrame = get(h,'JavaFrame');
% pause(0.1);
% set(jFrame,'Maximized',1);
% pause(0.1);
% warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

QF = Q1flash;
FRAME = 5*frame;
t0 = T(end)/FRAME;

for i = 1:FRAME
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    %contourf(xc,yc,QF(:,:,1,j),20);colormap(nclCM(484,100));
    %mesh(xc,yc,QF(:,:,1,j));axis([Xc(1),Xc(end),Yc(1),Yc(end),min(min(min(min(QF)))) - 0.5,max(max(max(max(QF)))) + 0.5]);colormap(cool);
    p = pcolor(xc,yc,QF(:,:,1,j));p.EdgeColor = 'none';colormap(jet);colorbar;
    %imshow(QF(:,:,j)'/max(max(max(abs(QF(:,:,j))))))
    %colormap(cool)
    pause(0.0001);
end
