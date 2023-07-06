    subroutine real_Ez
    
    include 'com.txt'
    
    include 'init1.txt'
    
    do i = 0,Nx
        do j = 1,Ny
            do j1 = 1,NumGLP
                Ez1real(i,j,j1) = U5(xa + i*hx - t,Yc(j) + hy1*lambda(j1) - t)
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do i1 = 1,NumGLP
                Ez2real(i,j,i1) = U6(Xc(i) + hx1*lambda(i1) - t,ya + j*hy - t)
            end do
        end do
    end do
    
    Bx = 0
    By = 0
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,k + 1
                do j1 = 1,NumGLP
                    Bx(i,j,d) = Bx(i,j,d) + 0.5*weight(j1)*EzG(j1,d)*Ez1real(i,j,j1)
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,k + 1
                do i1 = 1,NumGLP
                    By(i,j,d) = By(i,j,d) + 0.5*weight(i1)*EzG(i1,d)*Ez2real(i,j,i1)
                end do
            end do
        end do
    end do
    
    do d = 1,k + 1
        Bx(:,:,d) = Bx(:,:,d)/mmE(d)
        By(:,:,d) = By(:,:,d)/mmE(d)
    end do
    
    end subroutine real_Ez