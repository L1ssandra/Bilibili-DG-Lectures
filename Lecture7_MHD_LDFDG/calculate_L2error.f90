    subroutine calculate_L2error
    
    include 'com.txt'
    
    include 'init1.txt'
    
    L2 = 0
    uG = 0
    
    ! The value of num solution on the GL points
    do i = 1,Nx
        do j = 1,Ny
            uGint = 0
            do n = 1,NumEq
                do d = 1,dimPk
                    uGint(:,:,n) = uGint(:,:,n) + uh(i,j,d,n)*phiG(:,:,d)
                end do
            end do
            
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    L2(1) = L2(1) + weight(i1)*weight(j1)*(uGint(i1,j1,1) - U1(Xc(i) + hx1*lambda(i1) - t20,Yc(j) + hy1*lambda(j1) - t20))**2
                    L2(2) = L2(2) + weight(i1)*weight(j1)*(uGint(i1,j1,2) - U2(Xc(i) + hx1*lambda(i1) - t20,Yc(j) + hy1*lambda(j1) - t20))**2
                    L2(3) = L2(3) + weight(i1)*weight(j1)*(uGint(i1,j1,3) - U3(Xc(i) + hx1*lambda(i1) - t20,Yc(j) + hy1*lambda(j1) - t20))**2
                    L2(4) = L2(4) + weight(i1)*weight(j1)*(uGint(i1,j1,4) - U4(Xc(i) + hx1*lambda(i1) - t20,Yc(j) + hy1*lambda(j1) - t20))**2
                    L2(5) = L2(5) + weight(i1)*weight(j1)*(uGint(i1,j1,5) - U5(Xc(i) + hx1*lambda(i1) - t20,Yc(j) + hy1*lambda(j1) - t20))**2
                    L2(6) = L2(6) + weight(i1)*weight(j1)*(uGint(i1,j1,6) - U6(Xc(i) + hx1*lambda(i1) - t20,Yc(j) + hy1*lambda(j1) - t20))**2
                    L2(7) = L2(7) + weight(i1)*weight(j1)*(uGint(i1,j1,7) - U7(Xc(i) + hx1*lambda(i1) - t20,Yc(j) + hy1*lambda(j1) - t20))**2
                    L2(8) = L2(8) + weight(i1)*weight(j1)*(uGint(i1,j1,8) - U8(Xc(i) + hx1*lambda(i1) - t20,Yc(j) + hy1*lambda(j1) - t20))**2
                end do
            end do
        end do
    end do
    
    L2 = (L2*hx1*hy1)**0.5d0
    
    print *,L2
    
    end subroutine calculate_L2error