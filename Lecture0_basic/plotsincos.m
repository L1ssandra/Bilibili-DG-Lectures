% plotsincos.m

X1 = 0:0.01:2*pi;
Y1 = sin(X1);

X2 = 0:0.01:2*pi;
Y2 = cos(X2);

plot(X1,Y1,X2,Y2)
s = [0,2*pi,-1,1];
axis(s)