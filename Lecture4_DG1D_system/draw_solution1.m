%draw_solution1.m
figure(1)
hold on
plot(Xc,uh(:,1,1),'k-','linewidth',1.1)
plot(Xc1,u1(:,1,1),'r.-',Xc1,u2(:,1,1),'b.-',Xc1,u3(:,1,1),'c.-');
axis([0,1,min(uh(:,1)) - 0.2,max(uh(:,1)) + 0.2]);
legend('exact','LF','HLL','HLLC')

figure(2)
hold on
uhv1 = u1(:,1,2)./u1(:,1,1);
uhv2 = u2(:,1,2)./u2(:,1,1);
uhv3 = u3(:,1,2)./u3(:,1,1);
uhv = uh(:,1,2)./uh(:,1,1);
plot(Xc,uhv,'k-','linewidth',1.1)
plot(Xc1,uhv1,'r.-',Xc1,uhv2,'b.-',Xc1,uhv3,'c.-'); 
axis([0,1,min(uhv3) - 0.5,max(uhv3) + 0.5]);
legend('exact','LF','HLL','HLLC')

figure(3)
hold on
uhp1 = (gamma - 1)*(u1(:,1,3) - 0.5*u1(:,1,2).^2./u1(:,1,1));
uhp2 = (gamma - 1)*(u2(:,1,3) - 0.5*u2(:,1,2).^2./u2(:,1,1));
uhp3 = (gamma - 1)*(u3(:,1,3) - 0.5*u3(:,1,2).^2./u3(:,1,1));
uhp = (gamma - 1)*(uh(:,1,3) - 0.5*uh(:,1,2).^2./uh(:,1,1));
plot(Xc,uhp,'k-','linewidth',1.1)
plot(Xc1,uhp1,'r.-',Xc1,uhp2,'b.-',Xc1,uhp3,'c.-'); 
axis([0,1,min(uhp) - 50,max(uhp) + 50]);
legend('exact','LF','HLL','HLLC')