%draw_solution.m
figure(1)
plot(Xc,uh(:,1,1),'r.','linewidth',1.3); 
axis([0,1,min(uh(:,1,1)) - 0.5,max(uh(:,1,1)) + 0.5]);

figure(2)
uhv = uh(:,1,2)./uh(:,1,1);
plot(Xc,uhv,'r.','linewidth',1.3); 
axis([0,1,min(uhv) - 0.5,max(uhv) + 0.5]);

figure(3)
uhp = (gamma - 1)*(uh(:,1,3) - 0.5*uh(:,1,2).^2./uh(:,1,1));
plot(Xc,uhp,'r.','linewidth',1.3); 
axis([0,1,min(uhp) - 50,max(uhp) + 50]);