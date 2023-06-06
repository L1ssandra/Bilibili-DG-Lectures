    subroutine calculate_pmin
    
    include 'com.txt'
    
    pmin = 1000000
    
    do i = 1,Nx
        do j = 1,Ny
            uhGLL = 0
            do n = 1,NumEq
                do d = 1,dimPk
                    uhGLL(:,:,n,:) = uhGLL(:,:,n,:) + uh(i,j,d,n)*phiGLL(:,:,d,:)
                end do
            end do
            
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    do d = 1,3
                        p1 = pressure(uhGLL(i1,j1,1,d),uhGLL(i1,j1,2,d),uhGLL(i1,j1,3,d),uhGLL(i1,j1,4,d),gamma)
                        !p1 = pressure(uh(i,j,1,1),uh(i,j,1,2),uh(i,j,1,3),uh(i,j,1,4),uh(i,j,1,5),uh(i,j,1,6),uh(i,j,1,7),uh(i,j,1,8),gamma)
                        if (p1 < pmin) then
                            pmin = p1
                        end if
                    end do
                end do
            end do
            
        end do
    end do
    
    end subroutine calculate_pmin