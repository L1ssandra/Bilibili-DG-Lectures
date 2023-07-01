    subroutine init_data
    
    include 'com.txt'
    
    ! 0: Euler Vortex
    ! 1: MHD Vortex
    ! 2: Smooth Alfvén Wave 还没调好
    ! 3: Orszag-Tang Vortex
    ! 4: Rotor
    ! 5: MHD Blast
    ! 6: Cloud Shock Interaction
    ! 6b: Cloud Shock Interaction 2
    ! 7: Magnetic Loop
    ! 8: Brio-Wu Shock Tube
    ! 8b: Brio-Wu Shock Tube y-direction
    ! 9: Ryu-Jones Shock Tube
    ! 9b: Ryu-Jones Shock Tube y-direction
    
    include 'init4.txt'
    
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
                        uh(i,j,d,5) = uh(i,j,d,5) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U5(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                        uh(i,j,d,6) = uh(i,j,d,6) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U6(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                        uh(i,j,d,7) = uh(i,j,d,7) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U7(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                        uh(i,j,d,8) = uh(i,j,d,8) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U8(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                    end do
                end do
            end do
        end do
    end do
    
    uh = 0.25*uh
    
    do d = 1,dimPk
        uh(:,:,d,:) = uh(:,:,d,:)/mm(d)
    end do
    
    Bx = 0
    By = 0
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,k + 1
                do j1 = 1,NumGLP
                    Bx(i,j,d) = Bx(i,j,d) + 0.5*weight(j1)*EzG(j1,d)*U6(xa + i*hx,Yc(j) + hy1*lambda(j1))
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,k + 1
                do i1 = 1,NumGLP
                    By(i,j,d) = By(i,j,d) + 0.5*weight(i1)*EzG(i1,d)*U7(Xc(i) + hx1*lambda(i1),ya + j*hy)
                end do
            end do
        end do
    end do
    
    do d = 1,k + 1
        Bx(:,:,d) = Bx(:,:,d)/mmE(d)
        By(:,:,d) = By(:,:,d)/mmE(d)
    end do
    
    if (shock == 1) then
        do i = 1,Nx
            do j = 1,Ny
                uh(i,j,2:dimPk,:) = 0
                uh(i,j,1,1) = U1(Xc(i),Yc(j))
                uh(i,j,1,2) = U2(Xc(i),Yc(j))
                uh(i,j,1,3) = U3(Xc(i),Yc(j))
                uh(i,j,1,4) = U4(Xc(i),Yc(j))
                uh(i,j,1,5) = U5(Xc(i),Yc(j))
                uh(i,j,1,6) = U6(Xc(i),Yc(j))
                uh(i,j,1,7) = U7(Xc(i),Yc(j))
                uh(i,j,1,8) = U8(Xc(i),Yc(j))
            end do
        end do
        
        do i = 0,Nx
            do j = 1,Ny
                Bx(i,j,2:k + 1) = 0
                Bx(i,j,1) = U6(xa + i*hx,Yc(j))
            end do
        end do
        
        do i = 1,Nx
            do j = 0,Ny
                By(i,j,2:k + 1) = 0
                By(i,j,1) = U7(Xc(i),ya + j*hy)
            end do
        end do
        
    end if
    
    end subroutine init_data