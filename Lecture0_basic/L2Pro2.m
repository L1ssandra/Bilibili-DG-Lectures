% L2Pro2.m

f = @(x) sin(x);

p0 = @(x) 1 + 0.*x;
p1 = @(x) x;
p2 = @(x) x.^2 - 1/3;

xa = 0;
xb = 2*pi;
%N = 5;
h = (xb - xa)/N;
h1 = h/2;

m0 = 1;
m1 = 1/3;
m2 = 4/45;

lambda(1) = -0.8611363115940525752239465;
lambda(2) = -0.3399810435848562648026658;
lambda(3) = 0.3399810435848562648026658;
lambda(4) = 0.8611363115940525752239465;

weight(1) = 0.3478548451374538573730639;
weight(2) = 0.6521451548625461426269361;
weight(3) = 0.6521451548625461426269361;
weight(4) = 0.3478548451374538573730639;

Xc = xa + h1:h:xb - h1;
X = zeros(1,4*N);
for i = 1:N
    for j = 1:4
        X((i - 1)*4 + 1:i*4) = Xc(i) + h1*lambda;
    end
end
Y1 = f(X);
Y2 = zeros(1,4*N);

for i = 1:N
    I0 = 0;
    I1 = 0;
    I2 = 0;
    for j = 1:4
        I0 = I0 + weight(j)*f(Xc(i) + h1*lambda(j))*p0(lambda(j));
        I1 = I1 + weight(j)*f(Xc(i) + h1*lambda(j))*p1(lambda(j));
        I2 = I2 + weight(j)*f(Xc(i) + h1*lambda(j))*p2(lambda(j));
    end
    I0 = I0/(2*m0);
    I1 = I1/(2*m1);
    I2 = I2/(2*m2);
    fhi = @(x) I0*p0((x - Xc(i))/h1) + I1*p1((x - Xc(i))/h1) + I2*p2((x - Xc(i))/h1);
    for j = 1:4
        Y2((i - 1)*4 + j) = fhi(X((i - 1)*4 + j));
    end
end
%plot(X,Y1,'b-',X,Y2,'ro')
%s = [0,2*pi,-1.1,1.1];
%axis(s);
%legend('exact','num')

L2_Error = 0;
for i = 1:N
    for j = 1:4
        L2_Error = L2_Error + weight(j)*(Y2((i - 1)*4 + j) - Y1((i - 1)*4 + j))^2;
    end
end
L2_Error = sqrt(h1*L2_Error);