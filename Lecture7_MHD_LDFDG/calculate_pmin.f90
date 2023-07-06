    subroutine calculate_pmin
    
    include 'com.txt'
    
    pmin = 1000000
    
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
                    p1 = pressure(uGint(i1,j1,1),uGint(i1,j1,2),uGint(i1,j1,3),uGint(i1,j1,4),uGint(i1,j1,5),uGint(i1,j1,6),uGint(i1,j1,7),uGint(i1,j1,8),gamma)
                    if (p1 < pmin) then
                        pmin = p1
                    end if
                end do
            end do
            
        end do
    end do
    
    end subroutine calculate_pmin