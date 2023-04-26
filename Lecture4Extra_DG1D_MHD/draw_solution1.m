%draw_solution1.m
figure(1)
hold on
plot(Xc1,u1(:,1,1),'r.-',Xc1,u2(:,1,1),'m.-',Xc1,u3(:,1,1),'b.-',Xc1,u4(:,1,1),'c.-');
plot(Xc,uh(:,1,1),'k-','linewidth',1.1)
axis([-1,1,min(uh(:,1,1)) - 0.1,max(uh(:,1,1)) + 0.1]);
legend('LF','HLL','HLLC','HLLD','exact')

figure(2)
hold on
uhu1 = u1(:,1,2)./u1(:,1,1);
uhu2 = u2(:,1,2)./u2(:,1,1);
uhu3 = u3(:,1,2)./u3(:,1,1);
uhu4 = u4(:,1,2)./u4(:,1,1);
uhu = uh(:,1,2)./uh(:,1,1);
plot(Xc1,uhu1,'r.-',Xc1,uhu2,'m.-',Xc1,uhu3,'b.-',Xc1,uhu4,'c.-'); 
plot(Xc,uhu,'k-','linewidth',1.1)
axis([-1,1,min(uhu) - 0.1,max(uhu) + 0.1]);
legend('LF','HLL','HLLC','HLLD','exact')

figure(3)
hold on
uhv1 = u1(:,1,3)./u1(:,1,1);
uhv2 = u2(:,1,3)./u2(:,1,1);
uhv3 = u3(:,1,3)./u3(:,1,1);
uhv4 = u4(:,1,3)./u4(:,1,1);
uhv = uh(:,1,3)./uh(:,1,1);
plot(Xc1,uhv1,'r.-',Xc1,uhv2,'m.-',Xc1,uhv3,'b.-',Xc1,uhv4,'c.-'); 
plot(Xc,uhv,'k-','linewidth',1.1)
axis([-1,1,min(uhv) - 0.1,max(uhv) + 0.1]);
legend('LF','HLL','HLLC','HLLD','exact')

figure(4)
hold on
uhp1 = (gamma - 1)*(u1(:,1,4) - 0.5*(u1(:,1,2).^2 + u1(:,1,3).^2)./u1(:,1,1) - 0.5*(u1(:,1,5).^2 + u1(:,1,6).^2));
uhp2 = (gamma - 1)*(u2(:,1,4) - 0.5*(u2(:,1,2).^2 + u2(:,1,3).^2)./u2(:,1,1) - 0.5*(u2(:,1,5).^2 + u2(:,1,6).^2));
uhp3 = (gamma - 1)*(u3(:,1,4) - 0.5*(u3(:,1,2).^2 + u3(:,1,3).^2)./u3(:,1,1) - 0.5*(u3(:,1,5).^2 + u3(:,1,6).^2));
uhp4 = (gamma - 1)*(u4(:,1,4) - 0.5*(u4(:,1,2).^2 + u4(:,1,3).^2)./u4(:,1,1) - 0.5*(u4(:,1,5).^2 + u4(:,1,6).^2));
uhp = (gamma - 1)*(uh(:,1,4) - 0.5*(uh(:,1,2).^2 + uh(:,1,3).^2)./uh(:,1,1) - 0.5*(uh(:,1,5).^2 + uh(:,1,6).^2));
plot(Xc1,uhp1,'r.-',Xc1,uhp2,'m.-',Xc1,uhp3,'b.-',Xc1,uhp4,'c.-'); 
plot(Xc,uhp,'k-','linewidth',1.1)
axis([-1,1,min(uhp) - 0.1,max(uhp) + 0.1]);
legend('LF','HLL','HLLC','HLLD','exact')

figure(5)
hold on
plot(Xc1,u1(:,1,6),'r.-',Xc1,u2(:,1,6),'m.-',Xc1,u3(:,1,6),'b.-',Xc1,u4(:,1,6),'c.-');
plot(Xc,uh(:,1,6),'k-','linewidth',1.1)
axis([-1,1,min(uh(:,1,6)) - 0.1,max(uh(:,1,6)) + 0.1]);
legend('LF','HLL','HLLC','HLLD','exact')

