h = figure();				% ����ͼ�δ���
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	% �ر���صľ�����ʾ����Ϊ�����˷ǹ����ӿڣ�
jFrame = get(h,'JavaFrame');	% ��ȡ�ײ� Java �ṹ��ؾ��
pause(0.1);					% �� Win 10��Matlab 2017b �����²���ͣ�ٻᱨ Java �ײ���󡣸��˸�����Ҫ���Խ���ʵ����֤
set(jFrame,'Maximized',1);	%���������Ϊ�棨0 Ϊ�٣�
pause(0.1);					% ����ʵ���з��������ͣ�٣����ڿ����������仯������ȡ�Ĵ��ڴ�С����ԭ���ĳߴ硣���˸�����Ҫ���Խ���ʵ����֤
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		% ����ؾ�������

Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
QE = Q4;
QB1 = Q5;
QB2 = Q6;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP = 0.5*(QB1.^2 + QB2.^2);

% for i = 1:Nx
%     for j = 1:Ny
%         XY(i,j,:) = Xc(i) - Yc(j);
%     end
% end
% 
% Q1 = Q1 + XY;

TT = 30;
t0 = T(end)/TT;
%for k = 1:10
Qflash = QMach;

for i = 1:TT + 1
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    %mesh(yc,xc,Qflash(:,:,j));axis([Xc(1),Xc(end),Yc(1),Yc(end),min(min(min(Qflash))) - 0.01,max(max(max(Qflash))) + 0.01]);
    %colormap(cool);contour(yc,xc,Q1(:,:,j),15);
    p = pcolor(yc,xc,Qflash(:,:,j));colormap(jet);p.EdgeColor = 'none';
    title(T(j))
    pause(0.0001);
end
%end
% [xc,yc] = meshgrid(Xc(STOPL:STOPR),Yc(STOPL:STOPR));
% [~,j] = min(abs(T - 4));
% 
% QBPplot = QBP;
% Q1plot = Q1;

% for i = 1:N
%     for j = 1:N
%         if Q1plot(i,j,end) < 0.730 || Q1plot(i,j,end) > 7.330
%             Q1plot(i,j,end) = 0/0;
%         end
%         if QBPplot(i,j,end) < 0.059 || QBPplot(i,j,end) > 0.655
%             QBPplot(i,j,end) = 0/0;
%         end
%     end
% end
% 
% figure(2)
% colormap(hot)
% contourf(yc,xc,QBPplot(:,:,end),15);
% 
% figure(3)
% colormap(hot)
% contourf(yc,xc,Q1plot(:,:,end),15);