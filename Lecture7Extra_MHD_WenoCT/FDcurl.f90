    subroutine FDcurl(A,B,hx,hy,Nx,Ny)
    
    real A(Nx,Ny),B(Nx,Ny,2),AR2(Nx,Ny),AR1(Nx,Ny),AL1(Nx,Ny),AL2(Nx,Ny),hx,hy,dz
    integer:: Nx,Ny
    integer:: boundconditionX = 1
    integer:: boundconditionY = 3
    
    if (boundconditionX == 1) then
        
        AR2(Nx - 1:Nx,:) = A(1:2,:)
        AR2(1:Nx - 2,:) = A(3:Nx,:)
        
        AR1(Nx,:) = A(1,:)
        AR1(1:Nx - 1,:) = A(2:Nx,:)
        
        AL1(1,:) = A(Nx,:)
        AL1(2:Nx,:) = A(1:Nx - 1,:)
        
        AL2(1:2,:) = A(Nx - 1:Nx,:)
        AL2(3:Nx,:) = A(1:Nx - 2,:)
        
    else if (boundconditionX == 3) then
        
        dz = A(1,2) - A(1,1)
        
        AR2(:,Ny - 1:Ny) = A(:,Ny - 1:Ny) + 2*dz
        AR2(:,1:Ny - 2) = A(:,3:Ny)
        
        AR1(:,Ny) = A(:,Ny) + dz
        AR1(:,1:Ny - 1) = A(:,2:Ny)
        
        AL1(:,1) = A(:,1) - dz
        AL1(:,2:Ny) = A(:,1:Ny - 1)
        
        AL2(:,1:2) = A(:,1:2) - 2*dz
        AL2(:,3:Ny) = A(:,1:Ny - 2)
        
    else if (boundconditionX == 4) then
        
        dz = -hx
        
        AR2(Nx - 1:Nx,:) = A(1:2,:) + Nx*dz
        AR2(1:Nx - 2,:) = A(3:Ny,:)
        
        AR1(Nx,:) = A(1,:) + Nx*dz
        AR1(:,1:Ny - 1) = A(:,2:Ny)
        
        AL1(1,:) = A(Nx,:) - Nx*dz
        AL1(:,2:Ny) = A(:,1:Ny - 1)
        
        AL2(1:2,:) = A(Nx - 1:Nx,:) - Nx*dz
        AL2(:,3:Ny) = A(:,1:Ny - 2)
        
    end if
    
    B(:,:,2) = -(AL2 - 8*AL1 + 8*AR1 - AR2)/(12d0*hx)
    
    if (boundconditionY == 1) then
        
        AR2(:,Ny - 1:Ny) = A(:,1:2)
        AR2(:,1:Ny - 2) = A(:,3:Ny)
        
        AR1(:,Ny) = A(:,1)
        AR1(:,1:Ny - 1) = A(:,2:Ny)
        
        AL1(:,1) = A(:,Ny)
        AL1(:,2:Ny) = A(:,1:Ny - 1)
        
        AL2(:,1:2) = A(:,Ny - 1:Ny)
        AL2(:,3:Ny) = A(:,1:Ny - 2)
        
    else if (boundconditionY == 3) then
        
        dz = A(1,2) - A(1,1)
        
        AR2(:,Ny - 1:Ny) = A(:,Ny - 1:Ny) + 2*dz
        AR2(:,1:Ny - 2) = A(:,3:Ny)
        
        AR1(:,Ny) = A(:,Ny) + dz
        AR1(:,1:Ny - 1) = A(:,2:Ny)
        
        AL1(:,1) = A(:,1) - dz
        AL1(:,2:Ny) = A(:,1:Ny - 1)
        
        AL2(:,1:2) = A(:,1:2) - 2*dz
        AL2(:,3:Ny) = A(:,1:Ny - 2)
        
    else if (boundconditionY == 4) then
        
        dz = hy
        
        AR2(:,Ny - 1:Ny) = A(:,1:2) + Ny*dz
        AR2(:,1:Ny - 2) = A(:,3:Ny)
        
        AR1(:,Ny) = A(:,1) + Ny*dz
        AR1(:,1:Ny - 1) = A(:,2:Ny)
        
        AL1(:,1) = A(:,Ny) - Ny*dz
        AL1(:,2:Ny) = A(:,1:Ny - 1)
        
        AL2(:,1:2) = A(:,Ny - 1:Ny) - Ny*dz
        AL2(:,3:Ny) = A(:,1:Ny - 2)
        
    end if
    
    B(:,:,1) = (AL2 - 8*AL1 + 8*AR1 - AR2)/(12d0*hy)
    
    end subroutine FDcurl
    