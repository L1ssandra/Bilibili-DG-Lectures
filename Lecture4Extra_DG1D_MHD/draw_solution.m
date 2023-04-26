%draw_solution.m
figure(1)
plot(Xc,uh(:,1,1),'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uh(:,1,1)) - 0.1,max(uh(:,1,1)) + 0.1]);

figure(2)
uhu = uh(:,1,2)./uh(:,1,1);
plot(Xc,uhu,'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uhu) - 0.1,max(uhu) + 0.1]);

figure(3)
uhv = uh(:,1,3)./uh(:,1,1);
plot(Xc,uhv,'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uhv) - 0.1,max(uhv) + 0.1]);

figure(4)
uhp = (gamma - 1)*(uh(:,1,4) - 0.5*(uh(:,1,2).^2 + uh(:,1,3).^2)./uh(:,1,1) - 0.5*(uh(:,1,5).^2 + uh(:,1,6).^2));
plot(Xc,uhp,'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uhp) - 0.1,max(uhp) + 0.1]);

figure(5)
plot(Xc,uh(:,1,6),'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uh(:,1,6)) - 0.1,max(uh(:,1,6)) + 0.1]);