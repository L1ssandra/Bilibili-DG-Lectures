    subroutine LhEz
    
    include 'com.txt'
    
    ! calculate Ez on faces
    do i = 0,Nx
        do j = 1,Ny
            do j1 = 1,NumGLP
                SRmax = UR(i,j,j1,2)/UR(i,j,j1,1)
                SLmax = UL(i + 1,j,j1,2)/UL(i + 1,j,j1,1)
                ! L-F Flux
                SR = max(abs(SRmax),abs(SLmax))
               EzRL(i,j,j1) = -0.5*(FR(i,j,j1,7) + FL(i + 1,j,j1,7) - SR*(UL(i + 1,j,j1,7) - UR(i,j,j1,7)))
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do i1 = 1,NumGLP
                SRmax = UU(i,j,i1,3)/UU(i,j,i1,1)
                SLmax = UD(i,j + 1,i1,3)/UD(i,j + 1,i1,1)
                ! L-F Flux
                SR = max(abs(SRmax),abs(SLmax))
                EzUD(i,j,i1) = 0.5*(FU(i,j,i1,6) + FD(i,j + 1,i1,6) - SR*(UD(i,j + 1,i1,6) - UU(i,j,i1,6)))
            end do
        end do
    end do
    
    !EzRL = -Fxhat(:,:,:,7)
    !EzUD = Fyhat(:,:,:,6)
    
    ! calculate each Uh at vertex
    URU = 0
    ULU = 0
    URD = 0
    ULD = 0
    do i = 0,Nx1
        do j = 0,Ny1
            do d = 1,dimPk
                URU(i,j,:) = URU(i,j,:) + uh(i,j,d,:)*phiRU(d)
                ULU(i,j,:) = ULU(i,j,:) + uh(i,j,d,:)*phiLU(d)
                URD(i,j,:) = URD(i,j,:) + uh(i,j,d,:)*phiRD(d)
                ULD(i,j,:) = ULD(i,j,:) + uh(i,j,d,:)*phiLD(d)
            end do
        end do
    end do
    
    ! calculate Ez at vertex
    do i = 0,Nx
        do j = 0,Ny
            URU1 = ULD(i + 1,j + 1,:)
            ULU1 = URD(i,j + 1,:)
            URD1 = ULU(i + 1,j,:)
            ULD1 = URU(i,j,:)
            
            if (flux_type == 1) then
                call LF_Flux_2D
            else if (flux_type == 2) then
                call HLL_Flux_2D
            else if (flux_type == 3) then
                call HLLC_Flux_2D
            else if (flux_type == 4) then
                call HLLD_Flux_2D
            end if
            
            EzVertex(i,j) = Ezhat
        end do
    end do
    
    dEz1 = 0
    dEz2 = 0
    
    ! DG scheme of Ez1
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,k + 1
                do j1 = 1,NumGLP
                    dEz1(i,j,d) = dEz1(i,j,d) + 0.5*weight(j1)*EzRL(i,j,j1)*EzyG(j1,d)
                end do
                dEz1(i,j,d) = (dEz1(i,j,d) - EzVertex(i,j)*EzU(d)/hy + EzVertex(i,j - 1)*EzD(d)/hy)/mmE(d)
            end do
        end do
    end do
    
    ! DG scheme of Ez2
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,k + 1
                do i1 = 1,NumGLP
                    dEz2(i,j,d) = dEz2(i,j,d) - 0.5*weight(i1)*EzUD(i,j,i1)*EzxG(i1,d)
                end do
                dEz2(i,j,d) = (dEz2(i,j,d) + EzVertex(i,j)*EzR(d)/hx - EzVertex(i - 1,j)*EzL(d)/hx)/mmE(d)
            end do
        end do
    end do
    
    end subroutine LhEz