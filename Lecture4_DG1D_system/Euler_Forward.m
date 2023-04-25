% Euler_Forward.m

alpha = 1;
for i = 1:Nx
    [alpha1,~] = wavespeed(uh(i,:));
    if alpha1 > alpha
        alpha = alpha1;
    end
end

dt = CFL*hx/alpha;
t = 0;

while t < tend
    
    if t + dt >= tend
        dt = tend - t;
        t = tend;
    else
        t = t + dt;
    end
    
    du = Lh1(uh);
    uh = uh + dt*du;
    
    fprintf('%d  %d  %d\n',t,min(uh(:,1)),max(abs(uh(:,1))))
     
end