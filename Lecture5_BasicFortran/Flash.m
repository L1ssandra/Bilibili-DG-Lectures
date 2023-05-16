% Flash.m

% h = figure(1);
% warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
% jFrame = get(h,'JavaFrame');
% pause(0.1);
% set(jFrame,'Maximized',1);
% pause(0.1);
% warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

QF = uh;
FRAME = frame;
t0 = T(end)/(FRAME - 1);

for i = 1:FRAME
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    %contourf(x,y,QF(:,:,j),15);colormap(jet);
    %mesh(x,y,QF(:,:,j));axis([X(1),X(end),Y(1),Y(end),min(min(min(QF))),max(max(max(QF)))]);colormap(cool);
    p = pcolor(x,y,QF(:,:,j));colormap(jet);p.EdgeColor = 'none';
    colorbar;
    pause(0.0001);
end