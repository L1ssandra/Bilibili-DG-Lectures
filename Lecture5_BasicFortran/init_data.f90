    subroutine init_data
    
    include 'com.txt'
    
    include 'init1.txt'
    
    hx = (xb - xa)/Nx
    hy = (yb - ya)/Ny
    
    uh = 0
    
    open(unit = 1,file = 'X.txt')
    open(unit = 2,file = 'Y.txt')
    
    do i = 0,Nx
        X(i) = xa + i*hx
        write(1,*) X(i)
    end do
    
    do j = 0,Ny
        Y(j) = ya + j*hy
        write(2,*) Y(j)
    end do
    
    do i = 0,Nx
        do j = 0,Ny
            uh(i,j) = u0(X(i),Y(j))
        end do
    end do
    
    close(1)
    close(2)
    
    end subroutine init_data
    
    