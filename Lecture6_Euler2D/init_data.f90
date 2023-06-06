    subroutine init_data
    
    include 'com.txt'
    
    ! 0: Euler Vortex
    ! 1: sin
    ! 2: Sedov blast
    ! 3: Riemann problem
    ! 4: High Mach number jet
    ! 5: Double Mach reflection
    
    include 'init5.txt'
    
    hx = (xb - xa)/Nx
    hy = (yb - ya)/Ny
    hx1 = 0.5d0*hx
    hy1 = 0.5d0*hy
    
    do i = 1,Nx
        Xc(i) = xa + (i - 0.5)*hx
    end do
    
    do j = 1,Ny
        Yc(j) = ya + (j - 0.5)*hy
    end do
    
    call get_basis
    
    ! L2 Pro
    do i = 1,Nx
        do j = 1,Ny
            uh(i,j,:,:) = 0
            do d = 1,dimPk
                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        uh(i,j,d,1) = uh(i,j,d,1) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U1(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                        uh(i,j,d,2) = uh(i,j,d,2) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U2(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                        uh(i,j,d,3) = uh(i,j,d,3) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U3(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                        uh(i,j,d,4) = uh(i,j,d,4) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U4(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                    end do
                end do
            end do
        end do
    end do
    
    uh = 0.25*uh
    
    do d = 1,dimPk
        uh(:,:,d,:) = uh(:,:,d,:)/mm(d)
    end do
    
    if (shock == 1) then
        do i = 1,Nx
            do j = 1,Ny
                uh(i,j,2:dimPk,:) = 0
                uh(i,j,1,1) = U1(Xc(i),Yc(j))
                uh(i,j,1,2) = U2(Xc(i),Yc(j))
                uh(i,j,1,3) = U3(Xc(i),Yc(j))
                uh(i,j,1,4) = U4(Xc(i),Yc(j))
            end do
        end do
        
    end if
    
    if (blast == 1) then
        uh(1,1,1,4) = 0.244816/(hx*hy)
    end if
    
    end subroutine init_data