    subroutine pp_Limiter
    
    include 'com.txt'
    
    real eta,epsilon,eta1,pbar,pq,sq(NumEq),norm1,norm2,tq
    
    epsilon = 1e-13
    
    ! Limiting the density
    do i = 1,Nx
        do j = 1,Ny
            eta = 1
            uhGLL = 0
            do d = 1,dimPk
                do n = 1,NumEq
                    uhGLL(:,:,n,:) = uhGLL(:,:,n,:) + uh(i,j,d,n)*phiGLL(:,:,d,:)
                end do
            end do
            
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    do d = 1,2
                        ! eta1 = (rhobar - epsilon)/(rhobar - rho(xq))
                        eta1 = abs((uh(i,j,1,1) - epsilon)/(uh(i,j,1,1) - uhGLL(i1,j1,1,d)))
                        if (eta1 < 1) then
                            eta = eta1
                        end if
                    end do
                end do
            end do
            
            if (eta < 1) then
                uh(i,j,2:dimPk,1) = 0.5*eta*uh(i,j,2:dimPk,1)
            end if
            
        end do
    end do
    
    ! Limiting the pressure
    do i = 1,Nx
        do j = 1,Ny
            eta = 1
            eta1 = 1
            uhGLL = 0
            do d = 1,dimPk
                do n = 1,NumEq
                    uhGLL(:,:,n,:) = uhGLL(:,:,n,:) + uh(i,j,d,n)*phiGLL(:,:,d,:)
                end do
            end do
            pbar = pressure(uh(i,j,1,1),uh(i,j,1,2),uh(i,j,1,3),uh(i,j,1,4),uh(i,j,1,5),uh(i,j,1,6),uh(i,j,1,7),uh(i,j,1,8),gamma)
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    do d = 1,2
                        pq = pressure(uhGLL(i1,j1,1,d),uhGLL(i1,j1,2,d),uhGLL(i1,j1,3,d),uhGLL(i1,j1,4,d),uhGLL(i1,j1,5,d),uhGLL(i1,j1,6,d),uhGLL(i1,j1,7,d),uhGLL(i1,j1,8,d),gamma)
                        
                        if (pq < 0) then
                            call calculate_tq(uh(i,j,1,:),uhGLL(i1,j1,:,d),tq,gamma)
                            sq = tq*uhGLL(i1,j1,:,d) + (1 - tq)*uh(i,j,1,:)
                            call norm(sq - uh(i,j,1,:),norm1)
                            call norm(uhGLL(i1,j1,:,d) - uh(i,j,1,:),norm2)
                            eta1 = norm1/norm2
                            
                            !eta1 = pbar/(pbar - pq)
                        end if
                        
                        if (eta1 < eta) then
                            eta = eta1
                        end if
                        
                    end do
                end do
            end do
            
            if (eta < 1) then
                eta = 0.5*eta
            end if
            
            uh(i,j,2:dimPk,:) = eta*uh(i,j,2:dimPk,:)
            
        end do
    end do
    
    end subroutine pp_Limiter
    
    
    subroutine norm(x,d)
    
    real x(8),d
    
    d = 0
    
    do i = 1,8
        d = d + x(i)**2
    end do
    
    d = d**0.5
    
    end subroutine norm
    
    
    subroutine calculate_tq(ubar,uq,tq,gamma)
    
    real ubar(8),uq(8),tq,ta,tb,ut(8),gamma
    integer count
    
    ta = 0
    tb = 1
    count = 0
    
    do while (tb - ta > 1e-14)
        tq = 0.5*(ta + tb)
        ut = tq*uq + (1 - tq)*ubar
        if (pressure(ut(1),ut(2),ut(3),ut(4),ut(5),ut(6),ut(7),ut(8),gamma) < 1e-12) then
            !ta = ta
            tb = tq
        else
            ta = tq
            !tb = tb
        end if
    end do
    
    tq = ta
    ut = tq*uq + (1 - tq)*ubar
    !print *,pressure(ut(1),ut(2),ut(3),ut(4),ut(5),ut(6),ut(7),ut(8),gamma),ta,tb
    
    end subroutine calculate_tq