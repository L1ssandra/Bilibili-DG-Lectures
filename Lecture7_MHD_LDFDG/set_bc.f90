    subroutine set_bc
    
    include 'com.txt'
    
    if (bcR == 1) then
        uh(Nx1,:,:,:) = uh(1,:,:,:)
    else if (bcR == 2) then
        uh(Nx1,:,:,:) = uh(Nx,:,:,:)
    end if
    
    if (bcL == 1) then
        uh(0,:,:,:) = uh(Nx,:,:,:)
    else if (bcL == 2) then
        uh(0,:,:,:) = uh(1,:,:,:)
    end if
    
    if (bcU == 1) then
        uh(:,Ny1,:,:) = uh(:,1,:,:)
    else if (bcU == 2) then
        uh(:,Ny1,:,:) = uh(:,Ny,:,:)
    end if
    
    if (bcD == 1) then
        uh(:,0,:,:) = uh(:,Ny,:,:)
    else if (bcD == 2) then
        uh(:,0,:,:) = uh(:,1,:,:)
    end if
    
    end subroutine set_bc