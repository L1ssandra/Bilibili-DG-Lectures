    subroutine set_bc
    
    include 'com.txt'
    
    if (bcR == 1) then
        uh(Nx1,:,:,:) = uh(1,:,:,:)
    else if (bcR == 2) then
        uh(Nx1,:,:,:) = uh(Nx,:,:,:)
    end if
    
    if (bcL == 1) then
        uh(0,:,:,:) = uh(Nx,:,:,:)
    else if (bcL == 2) then
        uh(0,:,:,:) = uh(1,:,:,:)
    else if (bcL == 3) then
        uh(0,:,:,:) = 0
        uh(0,:,1,1) = 5
        do j = 1,Ny
            if (ya + (j - 0.5)*hy > -0.05 .and. ya + (j - 0.5)*hy < 0.05) then
                uh(0,j,1,2) = 5*800
                uh(0,j,1,4) = 0.4127/gamma1 + 0.5*5*800**2
            else
                uh(0,j,1,2) = 0
                uh(0,j,1,4) = 0.4127/gamma1
            end if
        end do
        uh(0,:,1,3) = 0
    else if (bcL == 4) then
        do j = 1,Ny
            call evenex_x(uh(0,j,:,1),uh(1,j,:,1))
            call oddex_x(uh(0,j,:,2),uh(1,j,:,2))
            call evenex_x(uh(0,j,:,3),uh(1,j,:,3))
            call evenex_x(uh(0,j,:,4),uh(1,j,:,4))
        end do
    end if
    
    if (bcU == 1) then
        uh(:,Ny1,:,:) = uh(:,1,:,:)
    else if (bcU == 2) then
        uh(:,Ny1,:,:) = uh(:,Ny,:,:)
    else if (bcU == 3) then
        uh(:,Ny1,:,:) = 0
        do i = 1,Nx
            if (Xc(i) < 1d0/6d0 + (1 + 20*tRK)/3d0**0.5) then ! post-shock
                uh(i,Ny1,1,1) = 8
                uh(i,Ny1,1,2) = 8*8.25*cos(pi/6d0)
                uh(i,Ny1,1,3) = -8*8.25*sin(pi/6d0)
                uh(i,Ny1,1,4) = 116.5d0/gamma1 + 0.5*8*((8.25*cos(pi/6d0))**2 + (8.25*sin(pi/6d0))**2)
            else ! pre-shock
                uh(i,Ny1,1,1) = 1.4
                uh(i,Ny1,1,2) = 0
                uh(i,Ny1,1,3) = 0
                uh(i,Ny1,1,4) = 1d0/gamma1
            end if
        end do
    end if
    
    if (bcD == 1) then
        uh(:,0,:,:) = uh(:,Ny,:,:)
    else if (bcD == 2) then
        uh(:,0,:,:) = uh(:,1,:,:)
    else if (bcD == 3) then
        do i = 1,Nx
            if (Xc(i) < 1d0/6d0) then
                uh(i,0,:,:) = uh(i,1,:,:)
            else
                call evenex_y(uh(i,0,:,1),uh(i,1,:,1))
                call evenex_y(uh(i,0,:,2),uh(i,1,:,2))
                call oddex_y(uh(i,0,:,3),uh(i,1,:,3))
                call evenex_y(uh(i,0,:,4),uh(i,1,:,4))
            end if
        end do
    else if (bcD == 4) then
        do i = 1,Nx
            call evenex_y(uh(i,0,:,1),uh(i,1,:,1))
            call evenex_y(uh(i,0,:,2),uh(i,1,:,2))
            call oddex_y(uh(i,0,:,3),uh(i,1,:,3))
            call evenex_y(uh(i,0,:,4),uh(i,1,:,4))
        end do
    end if
    
    end subroutine set_bc
    
    
    
    subroutine evenex_x(a,b)
    
    real a(6),b(6)
    
    a(1) = b(1)
    a(2) = -b(2)
    a(3) = b(3)
    a(4) = b(4)
    a(5) = -b(5)
    a(6) = b(6)
    
    end subroutine evenex_x
    
    
    subroutine oddex_x(a,b)
    
    real a(6),b(6)
    
    a(1) = -b(1)
    a(2) = b(2)
    a(3) = -b(3)
    a(4) = -b(4)
    a(5) = b(5)
    a(6) = -b(6)
    
    end subroutine oddex_x
    
    
    subroutine evenex_y(a,b)
    
    real a(6),b(6)
    
    a(1) = b(1)
    a(2) = b(2)
    a(3) = -b(3)
    a(4) = b(4)
    a(5) = -b(5)
    a(6) = b(6)
    
    end subroutine evenex_y
    
    
    subroutine oddex_y(a,b)
    
    real a(6),b(6)
    
    a(1) = -b(1)
    a(2) = -b(2)
    a(3) = b(3)
    a(4) = -b(4)
    a(5) = b(5)
    a(6) = -b(6)
    
    end subroutine oddex_y