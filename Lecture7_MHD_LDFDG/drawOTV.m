% drawOTV.m
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
%contour(rc,zc,Qrho,20);colormap(cool);
p = pcolor(rc,zc,Qrho);colormap(jet);p.EdgeColor = 'none';colorbar;
%mesh(xc,yc,Q1);
title('Density')

% figure(2)
% contour(xc,yc,QP,15);
% %mesh(xc,yc,Q1);
% %title('rho')
% colormap(cool);

% draw Bx at x = pi
figure(2)
[Nr,Nz] = size(Q6);
% Bx = load('Bx.txt');
% Bx = reshape(Bx,Nx/5,Ny);
% QBr = 0.5*(Q6(Nr/2,:) + Q6(Nr/2 + 1,:))./Rc;
plot(Zc,QBr,'k-','linewidth',1.3);
% plot(Yc,Bx(Nx/10,:),'r-','linewidth',1.3)
title('Br cut at r = 5\pi')
axis([Zc(1),Zc(end),min(QBr) - 0.2,max(QBr) + 0.2])