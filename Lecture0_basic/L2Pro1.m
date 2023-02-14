% L2Pro1.m

f = @(x) exp(x);

p0 = @(x) 1 + 0.*x;
p1 = @(x) x;
p2 = @(x) x.^2 - 1/3;

m0 = 2;
m1 = 2/3;
m2 = 8/45;

fp0 = @(x) f(x).*p0(x);
fp1 = @(x) f(x).*p1(x);
fp2 = @(x) f(x).*p2(x);

I0 = GaussIntgral(fp0,-1,1);
I1 = GaussIntgral(fp1,-1,1);
I2 = GaussIntgral(fp2,-1,1);

fh = @(x) (I0/m0)*p0(x) + (I1/m1)*p1(x) + (I2/m2)*p2(x);

X = -1:0.05:1;
Y1 = f(X);
Y2 = fh(X);
plot(X,Y1,'b-',X,Y2,'r-')
legend('exact','L2 Projection')