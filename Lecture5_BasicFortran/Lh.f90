    subroutine Lh
    
    include 'com.txt'
    
    du = 0
    
    call set_bc
    
    do i = 0,Nx
        do j = 0,Ny
            du(i,j) = -(uh(i,j) - uh(i - 1,j))/hx - (uh(i,j) - uh(i,j - 1))/hy
        end do
    end do
    
    end subroutine Lh