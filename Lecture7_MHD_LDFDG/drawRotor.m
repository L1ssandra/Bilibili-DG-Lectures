%drawRotor.m
Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
QE = Q5;
QB1 = Q6;
QB2 = Q7;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP = 0.5*(QB1.^2 + QB2.^2);

% for i = 1:Nx
%     for j = 1:Ny
%         if Qrho(i,j) < 1.127 || Qrho(i,j) > 5.587
%             Qrho(i,j) = NaN;
%         end
%     end
% end

figure(1);
contour(rc,zc,Q1,15);colormap(cool);
%p = pcolor(rc,zc,Q1);colormap(jet);p.EdgeColor = 'none';colorbar;
%mesh(xc,yc,Q1);
title('Density')

figure(2)
contour(rc,zc,QP,15);colormap(cool)
%p = pcolor(rc,zc,QP);colormap(jet);p.EdgeColor = 'none';colorbar;
%mesh(xc,yc,Q1);
%title('rho')
% colormap(cool);

figure(3)
contour(rc,zc,QMach,30);colormap(cool)
%p = pcolor(rc,zc,QMach);colormap(jet);p.EdgeColor = 'none';colorbar;
%mesh(xc,yc,Q1);
%title('rho')
% colormap(cool);

figure(4)
contour(rc,zc,QBP,15);colormap(cool)
%p = pcolor(rc,zc,QBP);colormap(jet);p.EdgeColor = 'none';colorbar;
%mesh(xc,yc,Q1);
%title('rho')
% colormap(cool);

% draw Br at x = 1.25
figure(5)
[Nr,Nz] = size(Q6);
% Bx = load('Bx.txt');
% Bx = reshape(Bx,Nx/5,Ny);
% QBr = 0.5*(Q6(Nr/2,:) + Q6(Nr/2 + 1,:))./Rc;
plot(Zc,QBrRotor,'k-','linewidth',1.3);
% plot(Yc,Bx(Nx/10,:),'r-','linewidth',1.3)
title('Br cut at r = 1.25')
axis([Zc(1),Zc(end),min(QBrRotor) - 0.2,max(QBrRotor) + 0.2])