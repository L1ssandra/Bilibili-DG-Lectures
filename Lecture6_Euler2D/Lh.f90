    subroutine Lh
    
    include 'com.txt'
    
    real rhoMR(0:Nx,Ny,NumGLP),uMR(0:Nx,Ny,NumGLP),vMR(0:Nx,Ny,NumGLP),wMR(0:Nx,Ny,NumGLP),EMR(0:Nx,Ny,NumGLP),pMR(0:Nx,Ny,NumGLP)
    real rhoML(Nx1,Ny,NumGLP),uML(Nx1,Ny,NumGLP),vML(Nx1,Ny,NumGLP),wML(Nx1,Ny,NumGLP),EML(Nx1,Ny,NumGLP),pML(Nx1,Ny,NumGLP)
    real rhoMU(Nx,0:Ny,NumGLP),uMU(Nx,0:Ny,NumGLP),vMU(Nx,0:Ny,NumGLP),wMU(Nx,0:Ny,NumGLP),EMU(Nx,0:Ny,NumGLP),pMU(Nx,0:Ny,NumGLP)
    real rhoMD(Nx,Ny1,NumGLP),uMD(Nx,Ny1,NumGLP),vMD(Nx,Ny1,NumGLP),wMD(Nx,Ny1,NumGLP),EMD(Nx,Ny1,NumGLP),pMD(Nx,Ny1,NumGLP)
    
    du = 0
    
    ! The value of num solution on the GL points
    
    !$omp parallel default(shared) private(uGint,rhoM,uM,vM,wM,EM,B1M,B2M,B3M,pM,SM,TM,KM,Fx,Fy,i,j,d,n,i1,j1)
    
    !$omp do collapse(2)
    do i = 1,Nx
        do j = 1,Ny
            
            uGint = 0
            do n = 1,NumEq
                do d = 1,dimPk
                    uGint(:,:,n) = uGint(:,:,n) + uh(i,j,d,n)*phiG(:,:,d)
                end do
            end do
            
            rhoM = uGint(:,:,1)
            uM = uGint(:,:,2)/rhoM
            vM = uGint(:,:,3)/rhoM
            EM = uGint(:,:,4)
    
            pM = gamma1*(EM - 0.5d0*rhoM*(uM**2 + vM**2))
    
            Fx(:,:,1) = uGint(:,:,2)
            Fx(:,:,2) = rhoM*uM**2 + pM
            Fx(:,:,3) = rhoM*uM*vM
            Fx(:,:,4) = uM*(EM + pM)
    
            Fy(:,:,1) = uGint(:,:,3)
            Fy(:,:,2) = rhoM*uM*vM
            Fy(:,:,3) = rhoM*vM**2 + pM
            Fy(:,:,4) = vM*(EM + pM)
            
            do d = 2,dimPk
                do n = 1,NumEq
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            du(i,j,d,n) = du(i,j,d,n) + 0.25d0*weight(i1)*weight(j1)*(Fx(i1,j1,n)*phixG(i1,j1,d) + Fy(i1,j1,n)*phiyG(i1,j1,d))
                        end do
                    end do
                end do
            end do
            
        end do
    end do
    !$omp end do
    
    ! The x-Flux
    !$omp single
    UR = 0
    UL = 0
    
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    UR(i,j,:,n) = UR(i,j,:,n) + uh(i,j,d,n)*phiGR(:,d)
                    UL(i + 1,j,:,n) = UL(i + 1,j,:,n) + uh(i + 1,j,d,n)*phiGL(:,d)
                end do
            end do
        end do
    end do
    !$omp end single
    
    !$omp workshare
    rhoMR = uR(:,:,:,1)
    uMR = uR(:,:,:,2)/rhoMR
    vMR = uR(:,:,:,3)/rhoMR
    EMR = uR(:,:,:,4)
    
    pMR = gamma1*(EMR - 0.5d0*rhoMR*(uMR**2 + vMR**2))
    
    FR(:,:,:,1) = uR(:,:,:,2)
    FR(:,:,:,2) = rhoMR*uMR**2 + pMR
    FR(:,:,:,3) = rhoMR*uMR*vMR
    FR(:,:,:,4) = uMR*(EMR + pMR)
    
    rhoML = uL(:,:,:,1)
    uML = uL(:,:,:,2)/rhoML
    vML = uL(:,:,:,3)/rhoML
    EML = uL(:,:,:,4)
    
    pML = gamma1*(EML - 0.5d0*rhoML*(uML**2 + vML**2))
    
    FL(:,:,:,1) = uL(:,:,:,2)
    FL(:,:,:,2) = rhoML*uML**2 + pML
    FL(:,:,:,3) = rhoML*uML*vML
    FL(:,:,:,4) = uML*(EML + pML)
    !$omp end workshare
    
    ! The y-Flux
    !$omp single
    UU = 0
    UD = 0
    
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    UU(i,j,:,n) = UU(i,j,:,n) + uh(i,j,d,n)*phiGU(:,d)
                    UD(i,j + 1,:,n) = UD(i,j + 1,:,n) + uh(i,j + 1,d,n)*phiGD(:,d)
                end do
            end do
        end do
    end do
    !$omp end single
    
    !$omp workshare
    rhoMU = UU(:,:,:,1)
    uMU = UU(:,:,:,2)/rhoMU
    vMU = UU(:,:,:,3)/rhoMU
    EMU = UU(:,:,:,4)
    
    pMU = gamma1*(EMU - 0.5d0*rhoMU*(uMU**2 + vMU**2))
    
    FU(:,:,:,1) = uU(:,:,:,3)
    FU(:,:,:,2) = rhoMU*uMU*vMU
    FU(:,:,:,3) = rhoMU*vMU**2 + pMU
    FU(:,:,:,4) = vMU*(EMU + pMU)
    
    rhoMD = UD(:,:,:,1)
    uMD = UD(:,:,:,2)/rhoMD
    vMD = UD(:,:,:,3)/rhoMD
    EMD = UD(:,:,:,4)
    
    pMD = gamma1*(EMD - 0.5d0*rhoMD*(uMD**2 + vMD**2))
    
    FD(:,:,:,1) = uD(:,:,:,3)
    FD(:,:,:,2) = rhoMD*uMD*vMD
    FD(:,:,:,3) = rhoMD*vMD**2 + pMD
    FD(:,:,:,4) = vMD*(EMD + pMD)
    !$omp end workshare
    
    ! calculate Fx hat
    !$omp single
    do i = 0,Nx
        do j = 1,Ny
            do j1 = 1,NumGLP
                call eigenvalueMm(SRmax,SRmin,UR(i,j,j1,1),UR(i,j,j1,2),UR(i,j,j1,3),UR(i,j,j1,4),1,0)
                call eigenvalueMm(SLmax,SLmin,UL(i + 1,j,j1,1),UL(i + 1,j,j1,2),UL(i + 1,j,j1,3),UL(i + 1,j,j1,4),1,0)
                !SR = 0.3*min(SRmax,SLmax)
                !SL = 0.3*max(SRmin,SLmin)
                SR = max(SRmax,SLmax)
                SL = min(SRmin,SLmin)
                FR1 = FL(i + 1,j,j1,:)
                FL1 = FR(i,j,j1,:)
                UR1 = UL(i + 1,j,j1,:)
                UL1 = UR(i,j,j1,:)
                if (flux_type == 1) then
                    call LF_Flux
                else if (flux_type == 2) then
                    call HLL_Flux
                else if (flux_type == 3) then
                    direction = 1
                    !call HLLC_Flux
                else if (flux_type == 4) then
                    direction = 1
                    !call HLLD_Flux
                end if
                Fxhat(i,j,j1,:) = Fhat1
            end do
        end do
    end do
    
    ! calculate Fy hat
    do i = 1,Nx
        do j = 0,Ny
            do i1 = 1,NumGLP
                call eigenvalueMm(SRmax,SRmin,UU(i,j,i1,1),UU(i,j,i1,2),UU(i,j,i1,3),UU(i,j,i1,4),0,1)
                call eigenvalueMm(SLmax,SLmin,UD(i,j + 1,i1,1),UD(i,j + 1,i1,2),UD(i,j + 1,i1,3),UD(i,j + 1,i1,4),0,1)
                !SR = 0.3*min(SRmax,SLmax)
                !SL = 0.3*max(SRmin,SLmin)
                SR = max(SRmax,SLmax)
                SL = min(SRmin,SLmin)
                FR1 = FD(i,j + 1,i1,:)
                FL1 = FU(i,j,i1,:)
                UR1 = UD(i,j + 1,i1,:)
                UL1 = UU(i,j,i1,:)
                if (flux_type == 1) then
                    call LF_Flux
                else if (flux_type == 2) then
                    call HLL_Flux
                else if (flux_type == 3) then
                    direction = 2
                    !call HLLC_Flux
                else if (flux_type == 4) then
                    direction = 2
                    !call HLLD_Flux
                end if
                Fyhat(i,j,i1,:) = Fhat1
            end do
        end do
    end do
    !$omp end single
    
    !!$omp single
    !$omp do collapse(5)
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    do j1 = 1,NumGLP
                        du(i,j,d,n) = du(i,j,d,n) - (0.5d0/hx)*weight(j1)*(Fxhat(i,j,j1,n)*phiGR(j1,d) - Fxhat(i - 1,j,j1,n)*phiGL(j1,d))
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    
    !$omp do collapse(5)
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    do i1 = 1,NumGLP
                        du(i,j,d,n) = du(i,j,d,n) - (0.5d0/hy)*weight(i1)*(Fyhat(i,j,i1,n)*phiGU(i1,d) - Fyhat(i,j - 1,i1,n)*phiGD(i1,d))
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    
    !$omp end parallel
    
    do d = 1,dimPk
        du(:,:,d,:) = du(:,:,d,:)/mm(d)
    end do
    
    end subroutine Lh