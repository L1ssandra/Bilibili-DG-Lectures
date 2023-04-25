% auto.m
global flux_type
flux_type = 1;
main
u1 = uh;
flux_type = 2;
main
u2 = uh;
flux_type = 3;
main
u3 = uh;
Xc1 = Xc;

uh1 = u1; uh2 = u2; uh3 = u3;

real_solution

figure(1)
hold on
plot(Xc1,uh1(:,1,1),'r.',Xc1,uh2(:,1,1),'b.',Xc1,uh3(:,1,1),'c.');
plot(Xc,uh(:,1),'k-','linewidth',1.1)
axis([0,1,min(uh(:,1)) - 0.5,max(uh(:,1)) + 0.5]);

figure(2)
hold on
uhv1 = uh1(:,1,2)./uh1(:,1,1);
uhv2 = uh2(:,1,2)./uh2(:,1,1);
uhv3 = uh3(:,1,2)./uh3(:,1,1);
uhv = uh(:,2)./uh(:,1);
plot(Xc1,uhv1,'r.',Xc1,uhv2,'b.',Xc1,uhv3,'c.'); 
plot(Xc,uhv,'k-','linewidth',1.1)
axis([0,1,min(uhv) - 0.5,max(uhv) + 0.5]);

figure(3)
hold on
uhp1 = (gamma - 1)*(uh1(:,1,3) - 0.5*uh1(:,1,2).^2./uh1(:,1,1));
uhp2 = (gamma - 1)*(uh2(:,1,3) - 0.5*uh2(:,1,2).^2./uh2(:,1,1));
uhp3 = (gamma - 1)*(uh3(:,1,3) - 0.5*uh3(:,1,2).^2./uh3(:,1,1));
uhp = (gamma - 1)*(uh(:,3) - 0.5*uh(:,2).^2./uh(:,1));
plot(Xc1,uhp1,'r.',Xc1,uhp2,'b.',Xc1,uhp3,'c.'); 
plot(Xc,uhp,'k-','linewidth',1.1)
axis([0,1,min(uhp) - 50,max(uhp) + 50]);