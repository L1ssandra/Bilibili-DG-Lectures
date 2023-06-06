% drawEz.m

Bx = load('Bx.txt');
By = load('By.txt');
Ez = load('Ez.txt');
Xc = load('Xc.txt');
Yc = load('Yc.txt');

Nx = length(Xc);
Ny = length(Yc);

xc = zeros(Nx,Ny);
yc = zeros(Nx,Ny);

for j = 1:Ny
    xc(:,j) = Xc;
end

for i = 1:Nx
    yc(i,:) = Yc';
end

Bx = reshape(Bx,Ny,Nx)';
By = reshape(By,Ny,Nx)';
Ez = reshape(Ez,Ny,Nx)';

figure(1)
mesh(xc,yc,Bx);colormap(bone)
figure(2)
mesh(xc,yc,By);colormap(bone)
figure(3)
mesh(xc,yc,Ez);colormap(bone)

