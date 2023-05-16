    subroutine set_bc
    
    include 'com.txt'
    
    if (bcR == 1) then
        uh(Nx1,:) = uh(1,:)
    end if
    
    if (bcL == 1) then
        uh(-1,:) = uh(Nx - 1,:)
    end if
    
    if (bcU == 1) then
        uh(:,Ny1) = uh(:,1)
    end if
    
    if (bcD == 1) then
        uh(:,-1) = uh(:,Ny - 1)
    end if
    
    end subroutine set_bc