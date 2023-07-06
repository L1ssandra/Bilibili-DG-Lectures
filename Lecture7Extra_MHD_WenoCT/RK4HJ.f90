subroutine RK4HJ(Q,QA,hx,hy,dt,Nx,Ny)
    
    ! 10-stages SSP-RK4
    integer i,dim
    parameter(dim = 1)
    real Q(Nx,Ny),Q1(Nx,Ny),Q2(Nx,Ny),hx,hy,dt,QLh(Nx,Ny),QA(Nx,Ny)
    
    Q1 = Q
    Q2 = Q
    do i = 1,5
        call Lh1(QLh,Q1,QA,hx,hy,Nx,Ny)
        Q1 = Q1 + (dt/6d0)*QLh
    end do
   
    Q2 = 0.04d0*Q2 + 0.36d0*Q1
    Q1 = 15*Q2 - 5*Q1
    do i = 6,9
        call Lh1(QLh,Q1,QA,hx,hy,Nx,Ny)
        Q1 = Q1 + (dt/6d0)*QLh
    end do
    
    call Lh1(QLh,Q1,QA,hx,hy,Nx,Ny)
    Q = Q2 + 0.6d0*Q1 + (dt/10d0)*QLh
    
    end subroutine RK4HJ
    
    
    subroutine Lh1(DQ,Q,QA,hx,hy,Nx,Ny)
    
    real DQ(Nx,Ny),Q(Nx,Ny),hx,hy,QR(Nx,Ny),QL(Nx,Ny),uR,uL,uU,uD,xu,yu,Hu,QA(Nx,Ny,6),dz
    real QR1(Nx,Ny),QR2(Nx,Ny),QR3(Nx,Ny),QL1(Nx,Ny),QL2(Nx,Ny),QL3(Nx,Ny)
    real QU1(Nx,Ny),QU2(Nx,Ny),QU3(Nx,Ny),QD1(Nx,Ny),QD2(Nx,Ny),QD3(Nx,Ny)
    integer Nx,Ny,i,j
    integer :: boundconditionX = 1
    integer :: boundconditionY = 3
    
    if (boundconditionX == 1) then
        
        QL3(1:3,:) = Q(Nx - 2:Nx,:)
        QL3(4:Nx,:) = Q(1:Nx - 3,:)
        
        QL2(1:2,:) = Q(Nx - 1:Nx,:)
        QL2(3:Nx,:) = Q(1:Nx - 2,:)
        
        QL1(1,:) = Q(Nx,:)
        QL1(2:Nx,:) = Q(1:Nx - 1,:)
        
        QR1(Nx,:) = Q(1,:)
        QR1(1:Nx - 1,:) = Q(2:Nx,:)
        
        QR2(Nx - 1:Nx,:) = Q(1:2,:)
        QR2(1:Nx - 2,:) = Q(3:Nx,:)
        
        QR3(Nx - 2:Nx,:) = Q(1:3,:)
        QR3(1:Nx - 3,:) = Q(4:Nx,:)
        
    else if (boundconditionX == 3) then
        
    else if (boundconditionX == 4) then
        
        dz = -hx
        
        QL3(1:3,:) = Q(Nx - 2:Nx,:) - Nx*dz
        QL3(4:Nx,:) = Q(1:Nx - 3,:)
        
        QL2(1:2,:) = Q(Nx - 1:Nx,:) - Nx*dz
        QL2(3:Nx,:) = Q(1:Nx - 2,:)
        
        QL1(1,:) = Q(Nx,:) - Nx*dz
        QL1(2:Nx,:) = Q(1:Nx - 1,:)
        
        QR1(Nx,:) = Q(1,:) + Nx*dz
        QR1(1:Nx - 1,:) = Q(2:Nx,:)
        
        QR2(Nx - 1:Nx,:) = Q(1:2,:) + Nx*dz
        QR2(1:Nx - 2,:) = Q(3:Nx,:)
        
        QR3(Nx - 2:Nx,:) = Q(1:3,:) + Nx*dz
        QR3(1:Nx - 3,:) = Q(4:Nx,:)
        
    end if
    
    if (boundconditionY == 1) then
        
        QD3(:,1:3) = Q(:,Ny - 2:Ny)
        QD3(:,4:Ny) = Q(:,1:Ny - 3)
        
        QD2(:,1:2) = Q(:,Ny - 1:Ny)
        QD2(:,3:Ny) = Q(:,1:Ny - 2)
        
        QD1(:,1) = Q(:,Ny)
        QD1(:,2:Ny) = Q(:,1:Ny - 1)
        
        QU1(:,Ny) = Q(:,1)
        QU1(:,1:Ny - 1) = Q(:,2:Ny)
        
        QU2(:,Ny - 1:Ny) = Q(:,1:2)
        QU2(:,1:Ny - 2) = Q(:,3:Ny)
        
        QU3(:,Ny - 2:Ny) = Q(:,1:3)
        QU3(:,1:Ny - 3) = Q(:,4:Ny)
    
    else if (boundconditionY == 3) then
        ! ÏßÐÔ±ß½ç
        dz = Q(1,2) - Q(1,1)
        
        QD3(:,1:3) = Q(:,1:3) - 3*dz
        QD3(:,4:Ny) = Q(:,1:Ny - 3)
        
        QD2(:,1:2) = Q(:,1:2) - 2*dz
        QD2(:,3:Ny) = Q(:,1:Ny - 2)
        
        QD1(:,1) = Q(:,1) - dz
        QD1(:,2:Ny) = Q(:,1:Ny - 1)
        
        QU1(:,Ny) = Q(:,Ny) + dz
        QU1(:,1:Ny - 1) = Q(:,2:Ny)
        
        QU2(:,Ny - 1:Ny) = Q(:,Ny - 1:Ny) + 2*dz
        QU2(:,1:Ny - 2) = Q(:,3:Ny)
        
        QU3(:,Ny - 2:Ny) = Q(:,Ny - 2:Ny) + 3*dz
        QU3(:,1:Ny - 3) = Q(:,4:Ny)
        
    else if (boundconditionY == 4) then
        
        dz = hy
        
        QD3(:,1:3) = Q(:,Ny - 2:Ny) - Ny*dz
        QD3(:,4:Ny) = Q(:,1:Ny - 3)
        
        QD2(:,1:2) = Q(:,Ny - 1:Ny) - Ny*dz
        QD2(:,3:Ny) = Q(:,1:Ny - 2)
        
        QD1(:,1) = Q(:,Ny) - Ny*dz
        QD1(:,2:Ny) = Q(:,1:Ny - 1)
        
        QU1(:,Ny) = Q(:,1) + Ny*dz
        QU1(:,1:Ny - 1) = Q(:,2:Ny)
        
        QU2(:,Ny - 1:Ny) = Q(:,1:2) + Ny*dz
        QU2(:,1:Ny - 2) = Q(:,3:Ny)
        
        QU3(:,Ny - 2:Ny) = Q(:,1:3) + Ny*dz
        QU3(:,1:Ny - 3) = Q(:,4:Ny)
        
    end if
    
    alphax = 0
    alphay = 0
    
    do i = 1,Nx
        do j = 1,Ny
            
            if (abs(QA(i,j,2)/QA(i,j,1)) > alphax) then
                alphax = abs(QA(i,j,2)/QA(i,j,1))
            end if
            
            if (abs(QA(i,j,3)/QA(i,j,1)) > alphay) then
                alphay = abs(QA(i,j,3)/QA(i,j,1))
            end if
            
        end do
    end do
            
    
    do i = 1,Nx
        do j = 1,Ny
            
            call WENO5( uL,(QL2(i,j) - QL3(i,j))/hx,(QL1(i,j) - QL2(i,j))/hx,(Q(i,j) - QL1(i,j))/hx,(QR1(i,j) - Q(i,j))/hx,(QR2(i,j) - QR1(i,j))/hx )
            call WENO5( uR,(QR3(i,j) - QR2(i,j))/hx,(QR2(i,j) - QR1(i,j))/hx,(QR1(i,j) - Q(i,j))/hx,(Q(i,j) - QL1(i,j))/hx,(QL1(i,j) - QL2(i,j))/hx )
            call WENO5( uD,(QD2(i,j) - QD3(i,j))/hy,(QD1(i,j) - QD2(i,j))/hy,(Q(i,j) - QD1(i,j))/hy,(QU1(i,j) - Q(i,j))/hy,(QU2(i,j) - QU1(i,j))/hy )
            call WENO5( uU,(QU3(i,j) - QU2(i,j))/hy,(QU2(i,j) - QU1(i,j))/hy,(QU1(i,j) - Q(i,j))/hy,(Q(i,j) - QD1(i,j))/hy,(QD1(i,j) - QD2(i,j))/hy )
            
            xu = 0.5d0*(uR + uL)
            yu = 0.5d0*(uU + uD)
            
            call H(Hu,xu,yu,QA(i,j,2)/QA(i,j,1),QA(i,j,3)/QA(i,j,1))
            
            DQ(i,j) = -(Hu - alphax*0.5d0*(uR - uL) - alphay*0.5d0*(uU - uD))
            
        end do
    end do
    
    end subroutine Lh1
    
    
    
    subroutine H(z,x,y,u1,u2)
    
    real x,y,u1,u2,z
    
    z = u1*x + u2*y
    
    end subroutine H
    