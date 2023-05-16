% readFortran.m

X = load('X.txt');
Y = load('Y.txt');
T = load('T.txt');
Nx = length(X) - 1;
Ny = length(Y) - 1;
Nx1 = Nx + 1;
Ny1 = Ny + 1;
frame = length(T);

uh = load('uh.txt');

x = zeros(Nx1,Ny1);
y = zeros(Nx1,Ny1);

for j = 1:Nx1
    x(:,j) = X;
end

for i = 1:Ny1
    y(i,:) = Y';
end

uh = reshape(uh,Nx1,Ny1,frame);