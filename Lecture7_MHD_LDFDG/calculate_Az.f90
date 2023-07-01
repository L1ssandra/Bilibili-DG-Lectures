    subroutine calculate_Az
    
    include 'com.txt'
    
    include 'init7.txt'
    
    do i = 0,Nx
        Xbb(i) = xa + i*hx
    end do
    
    do j = 0,Ny
        Ybb(j) = ya + j*hy
    end do
    
    do i = 0,Nx
        do j = 1,Ny
            Bx(i,j,1) = -(1d0/hy)*(Az(Xbb(i),Ybb(j)) - Az(Xbb(i),Ybb(j - 1)))
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            By(i,j,1) = (1d0/hx)*(Az(Xbb(i),Ybb(j)) - Az(Xbb(i - 1),Ybb(j)))
        end do
    end do
    
    end subroutine calculate_Az