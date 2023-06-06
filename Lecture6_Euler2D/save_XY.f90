    subroutine save_XY
    
    include 'com.txt'
    
    open(unit = 1,file = 'Xc.txt')
    open(unit = 2,file = 'Yc.txt')
    open(unit = 3,file = 'lambda.txt')
    open(unit = 4,file = 'weight.txt')
    
    do i = 1,Nx
        do i1 = 1,NumGLP
            write(1,*) Xc(i) + hx1*lambda(i1)
        end do
    end do
    
    do j = 1,Ny
        do j1 = 1,NumGLP
            write(2,*) Yc(j) + hy1*lambda(j1)
        end do
    end do
    
    do i = 1,NumGLP
        write(3,*) lambda(i)
        write(4,*) weight(i)
    end do
    
    close(1)
    close(2)
    close(3)
    close(4)
    
    end subroutine save_XY