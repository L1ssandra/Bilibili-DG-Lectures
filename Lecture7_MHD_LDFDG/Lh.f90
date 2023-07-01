    subroutine Lh
    
    include 'com.txt'
    
    du = 0
    
    ! The value of num solution on the GL points
    
    !$omp parallel num_threads(16) default(shared) private(uGint,rhoM,uM,vM,wM,EM,B1M,B2M,B3M,pM,SM,TM,KM,Fx,Fy,i,j,d,n,i1,j1)
    
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
            wM = uGint(:,:,4)/rhoM
            EM = uGint(:,:,5)
            B1M = uGint(:,:,6)
            B2M = uGint(:,:,7)
            B3M = uGint(:,:,8)
    
            pM = gamma1*(EM - 0.5d0*rhoM*(uM**2 + vM**2 + wM**2) - 0.5d0*(B1M**2 + B2M**2 + B3M**2))
    
            SM = pM + 0.5d0*(B1M**2 + B2M**2 + B3M**2)
            TM = EM + SM
            KM = uM*B1M + vM*B2M + wM*B3M
    
            Fx(:,:,1) = uGint(:,:,2)
            Fx(:,:,2) = rhoM*uM**2 + SM - B1M**2
            Fx(:,:,3) = rhoM*uM*vM - B1M*B2M
            Fx(:,:,4) = rhoM*uM*wM - B1M*B3M
            Fx(:,:,5) = TM*uM - KM*B1M
            Fx(:,:,6) = 0
            Fx(:,:,7) = uM*B2M - vM*B1M
            Fx(:,:,8) = uM*B3M - wM*B1M
    
            Fy(:,:,1) = uGint(:,:,3)
            Fy(:,:,2) = rhoM*uM*vM - B1M*B2M
            Fy(:,:,3) = rhoM*vM**2 + SM - B2M**2
            Fy(:,:,4) = rhoM*vM*wM - B2M*B3M
            Fy(:,:,5) = TM*vM - KM*B2M
            Fy(:,:,6) = vM*B1M - uM*B2M
            Fy(:,:,7) = 0
            Fy(:,:,8) = vM*B3M - wM*B2M
            
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
    wMR = uR(:,:,:,4)/rhoMR
    EMR = uR(:,:,:,5)
    B1MR = uR(:,:,:,6)
    B2MR = uR(:,:,:,7)
    B3MR = uR(:,:,:,8)
    
    pMR = gamma1*(EMR - 0.5d0*rhoMR*(uMR**2 + vMR**2 + wMR**2) - 0.5d0*(B1MR**2 + B2MR**2 + B3MR**2))
    
    SMR = pMR + 0.5d0*(B1MR**2 + B2MR**2 + B3MR**2)
    TMR = EMR + SMR
    KMR = uMR*B1MR + vMR*B2MR + wMR*B3MR
    
    FR(:,:,:,1) = uR(:,:,:,2)
    FR(:,:,:,2) = rhoMR*uMR**2 + SMR - B1MR**2
    FR(:,:,:,3) = rhoMR*uMR*vMR - B1MR*B2MR
    FR(:,:,:,4) = rhoMR*uMR*wMR - B1MR*B3MR
    FR(:,:,:,5) = TMR*uMR - KMR*B1MR
    FR(:,:,:,6) = 0
    FR(:,:,:,7) = uMR*B2MR - vMR*B1MR
    FR(:,:,:,8) = uMR*B3MR - wMR*B1MR
    
    rhoML = uL(:,:,:,1)
    uML = uL(:,:,:,2)/rhoML
    vML = uL(:,:,:,3)/rhoML
    wML = uL(:,:,:,4)/rhoML
    EML = uL(:,:,:,5)
    B1ML = uL(:,:,:,6)
    B2ML = uL(:,:,:,7)
    B3ML = uL(:,:,:,8)
    
    pML = gamma1*(EML - 0.5d0*rhoML*(uML**2 + vML**2 + wML**2) - 0.5d0*(B1ML**2 + B2ML**2 + B3ML**2))
    
    SML = pML + 0.5d0*(B1ML**2 + B2ML**2 + B3ML**2)
    TML = EML + SML
    KML = uML*B1ML + vML*B2ML + wML*B3ML
    
    FL(:,:,:,1) = uL(:,:,:,2)
    FL(:,:,:,2) = rhoML*uML**2 + SML - B1ML**2
    FL(:,:,:,3) = rhoML*uML*vML - B1ML*B2ML
    FL(:,:,:,4) = rhoML*uML*wML - B1ML*B3ML
    FL(:,:,:,5) = TML*uML - KML*B1ML
    FL(:,:,:,6) = 0
    FL(:,:,:,7) = uML*B2ML - vML*B1ML
    FL(:,:,:,8) = uML*B3ML - wML*B1ML
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
    wMU = UU(:,:,:,4)/rhoMU
    EMU = UU(:,:,:,5)
    B1MU = UU(:,:,:,6)
    B2MU = UU(:,:,:,7)
    B3MU = UU(:,:,:,8)
    
    pMU = gamma1*(EMU - 0.5d0*rhoMU*(uMU**2 + vMU**2 + wMU**2) - 0.5d0*(B1MU**2 + B2MU**2 + B3MU**2))
    
    SMU = pMU + 0.5d0*(B1MU**2 + B2MU**2 + B3MU**2)
    TMU = EMU + SMU
    KMU = uMU*B1MU + vMU*B2MU + wMU*B3MU
    
    FU(:,:,:,1) = uU(:,:,:,3)
    FU(:,:,:,2) = rhoMU*uMU*vMU - B1MU*B2MU
    FU(:,:,:,3) = rhoMU*vMU**2 + SMU - B2MU**2
    FU(:,:,:,4) = rhoMU*vMU*wMU - B2MU*B3MU
    FU(:,:,:,5) = TMU*vMU - KMU*B2MU
    FU(:,:,:,6) = vMU*B1MU - uMU*B2MU
    FU(:,:,:,7) = 0
    FU(:,:,:,8) = vMU*B3MU - wMU*B2MU
    
    rhoMD = UD(:,:,:,1)
    uMD = UD(:,:,:,2)/rhoMD
    vMD = UD(:,:,:,3)/rhoMD
    wMD = UD(:,:,:,4)/rhoMD
    EMD = UD(:,:,:,5)
    B1MD = UD(:,:,:,6)
    B2MD = UD(:,:,:,7)
    B3MD = UD(:,:,:,8)
    
    pMD = gamma1*(EMD - 0.5d0*rhoMD*(uMD**2 + vMD**2 - wMD**2) - 0.5d0*(B1MD**2 + B2MD**2 + B3MD**2))
    
    SMD = pMD + 0.5d0*(B1MD**2 + B2MD**2 + B3MD**2)
    TMD = EMD + SMD
    KMD = uMD*B1MD + vMD*B2MD + wMD*B3MD
    
    FD(:,:,:,1) = uD(:,:,:,3)
    FD(:,:,:,2) = rhoMD*uMD*vMD - B1MD*B2MD
    FD(:,:,:,3) = rhoMD*vMD**2 + SMD - B2MD**2
    FD(:,:,:,4) = rhoMD*vMD*wMD - B2MD*B3MD
    FD(:,:,:,5) = TMD*vMD - KMD*B2MD
    FD(:,:,:,6) = vMD*B1MD - uMD*B2MD
    FD(:,:,:,7) = 0
    FD(:,:,:,8) = vMD*B3MD - wMD*B2MD
    !$omp end workshare
    
    ! calculate Fx hat
    !$omp single
    do i = 0,Nx
        do j = 1,Ny
            do j1 = 1,NumGLP
                call eigenvalueMm(SRmax,SRmin,UR(i,j,j1,1),UR(i,j,j1,2),UR(i,j,j1,3),UR(i,j,j1,4),UR(i,j,j1,5),UR(i,j,j1,6),UR(i,j,j1,7),UR(i,j,j1,8),1,0)
                call eigenvalueMm(SLmax,SLmin,UL(i + 1,j,j1,1),UL(i + 1,j,j1,2),UL(i + 1,j,j1,3),UL(i + 1,j,j1,4),UL(i + 1,j,j1,5),UL(i + 1,j,j1,6),UL(i + 1,j,j1,7),UL(i + 1,j,j1,8),1,0)
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
                    call HLLC_Flux
                else if (flux_type == 4) then
                    direction = 1
                    call HLLD_Flux
                end if
                Fxhat(i,j,j1,:) = Fhat1
            end do
        end do
    end do
    
    ! calculate Fy hat
    do i = 1,Nx
        do j = 0,Ny
            do i1 = 1,NumGLP
                call eigenvalueMm(SRmax,SRmin,UU(i,j,i1,1),UU(i,j,i1,2),UU(i,j,i1,3),UU(i,j,i1,4),UU(i,j,i1,5),UU(i,j,i1,6),UU(i,j,i1,7),UU(i,j,i1,8),0,1)
                call eigenvalueMm(SLmax,SLmin,UD(i,j + 1,i1,1),UD(i,j + 1,i1,2),UD(i,j + 1,i1,3),UD(i,j + 1,i1,4),UD(i,j + 1,i1,5),UD(i,j + 1,i1,6),UD(i,j + 1,i1,7),UD(i,j + 1,i1,8),0,1)
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
                    call HLLC_Flux
                else if (flux_type == 4) then
                    direction = 2
                    call HLLD_Flux
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