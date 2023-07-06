    subroutine RK4Chara(Q,QA,hx,hy,dt,Nx,Ny)
    
    ! Ketcheson 10-stages SSP-RK4
    integer i,dim
    parameter(dim = 6)
    real Q(Nx,Ny,dim),Q1(Nx,Ny,dim),Q2(Nx,Ny,dim),hx,hy,dt,QLh(Nx,Ny,dim),ALh(Nx,Ny)
    real QA(Nx,Ny),B(Nx,Ny,2),QALh(Nx,Ny)
    real QA1(Nx,Ny),QA2(Nx,Ny)
    
    Q1 = Q
    Q2 = Q
    QA1 = QA
    QA2 = QA
    do i = 1,5
        call Lh(QLh,Q1,hx,hy,Nx,Ny)
        call Lh1(QALh,QA1,Q1,hx,hy,Nx,Ny)
        Q1 = Q1 + (dt/6d0)*QLh
        QA1 = QA1 + (dt/6d0)*QALh
        call FDcurl(QA1,B,hx,hy,Nx,Ny)
        Q1(:,:,5:6) = B
    end do
   
    Q2 = 0.04d0*Q2 + 0.36d0*Q1
    QA2 = 0.04d0*QA2 + 0.36d0*QA1
    Q1 = 15*Q2 - 5*Q1
    QA1 = 15*QA2 - 5*QA1
    do i = 6,9
        call Lh(QLh,Q1,hx,hy,Nx,Ny)
        call Lh1(QALh,QA1,Q1,hx,hy,Nx,Ny)
        Q1 = Q1 + (dt/6d0)*QLh
        QA1 = QA1 + (dt/6d0)*QALh
        call FDcurl(QA1,B,hx,hy,Nx,Ny)
        Q1(:,:,5:6) = B
    end do
    
    call Lh(QLh,Q1,hx,hy,Nx,Ny)
    call Lh1(QALh,QA1,Q1,hx,hy,Nx,Ny)
    Q = Q2 + 0.6d0*Q1 + (dt/10d0)*QLh
    QA = QA2 + 0.6*QA1 + (dt/10d0)*QALh
    call FDcurl(QA,B,hx,hy,Nx,Ny)
    Q(:,:,5:6) = B
    
    end subroutine RK4Chara
    
    
    subroutine Lh(DQ,Q,hx,hy,Nx,Ny)
    
    integer dim,dim1
    
    parameter(dim = 6)
    parameter(dim1 = 8)
    
    real Q(Nx,Ny,dim),hx,hy,DQ(Nx,Ny,dim),fQ(Nx,Ny,dim),yv(dim),DQx(Nx,Ny,dim),DQy(Nx,Ny,dim),fhatx(Nx + 1,Ny,dim),fhaty(Nx,Ny + 1,dim)
    real fQR1(Nx,Ny,dim),fQR2(Nx,Ny,dim),fQR3(Nx,Ny,dim),fQL1(Nx,Ny,dim),fQL2(Nx,Ny,dim),rho,u,v,p,E,H,VV,gamma1,B1,B2
    real QR1(Nx,Ny,dim),QR2(Nx,Ny,dim),QR3(Nx,Ny,dim),QL1(Nx,Ny,dim),QL2(Nx,Ny,dim),QRnew(dim1,1),QLnew(dim1,1)
    real fQP(dim1,5),fQN(dim1,5),alpha,AlphaV(dim1,1),AAlpha(dim1,5),alpha1,alpha2,alpha3,QP(dim1,5),QN(dim1,5),gamma,AlphaV1(dim1,1)
    real L(dim1,dim1),R(dim1,dim1)
    integer i,j
    integer :: boundcondition = 1
    integer :: undoeig = 0
    
    gamma = 5d0/3d0
    
    ! x-direction
    boundcondition = 1
    
    ! Step 1£ºcalculate f(q)
    do i = 1,Nx
        do j = 1,Ny
            call f1(Q(i,j,:),yv)
            fQ(i,j,:) = yv
        end do
    end do
    
    ! Step 2£ºset boundary condition
    
    if (boundcondition == 1) then
        
        ! f
        fQL2(1:2,:,:) = fQ(Nx - 1:Nx,:,:)
        fQL2(3:Nx,:,:) = fQ(1:Nx - 2,:,:)
        
        fQL1(1,:,:) = fQ(Nx,:,:)
        fQL1(2:Nx,:,:) = fQ(1:Nx - 1,:,:)
        
        fQR1(Nx,:,:) = fQ(1,:,:)
        fQR1(1:Nx - 1,:,:) = fQ(2:Nx,:,:)
        
        fQR2(Nx - 1:Nx,:,:) = fQ(1:2,:,:)
        fQR2(1:Nx - 2,:,:) = fQ(3:Nx,:,:)
        
        fQR3(Nx - 2:Nx,:,:) = fQ(1:3,:,:)
        fQR3(1:Nx - 3,:,:) = fQ(4:Nx,:,:)
        
        ! u
        QL2(1:2,:,:) = Q(Nx - 1:Nx,:,:)
        QL2(3:Nx,:,:) = Q(1:Nx - 2,:,:)
        
        QL1(1,:,:) = Q(Nx,:,:)
        QL1(2:Nx,:,:) = Q(1:Nx - 1,:,:)
        
        QR1(Nx,:,:) = Q(1,:,:)
        QR1(1:Nx - 1,:,:) = Q(2:Nx,:,:)
        
        QR2(Nx - 1:Nx,:,:) = Q(1:2,:,:)
        QR2(1:Nx - 2,:,:) = Q(3:Nx,:,:)
        
        QR3(Nx - 2:Nx,:,:) = Q(1:3,:,:)
        QR3(1:Nx - 3,:,:) = Q(4:Nx,:,:)
        
    else if (boundcondition == 2) then
        
        ! f
        fQL2(1:2,:,:) = fQ(1:2,:,:)
        fQL2(3:Nx,:,:) = fQ(1:Nx - 2,:,:)
        
        fQL1(1,:,:) = fQ(1,:,:)
        fQL1(2:Nx,:,:) = fQ(1:Nx - 1,:,:)
        
        fQR1(Nx,:,:) = fQ(Nx,:,:)
        fQR1(1:Nx - 1,:,:) = fQ(2:Nx,:,:)
        
        fQR2(Nx - 1:Nx,:,:) = fQ(Nx - 1:Nx,:,:)
        fQR2(1:Nx - 2,:,:) = fQ(3:Nx,:,:)
        
        fQR3(Nx - 2:Nx,:,:) = fQ(Nx - 2:Nx,:,:)
        fQR3(1:Nx - 3,:,:) = fQ(4:Nx,:,:)
        
        ! u
        QL2(1:2,:,:) = Q(1:2,:,:)
        QL2(3:Nx,:,:) = Q(1:Nx - 2,:,:)
        
        QL1(1,:,:) = Q(1,:,:)
        QL1(2:Nx,:,:) = Q(1:Nx - 1,:,:)
        
        QR1(Nx,:,:) = Q(Nx,:,:)
        QR1(1:Nx - 1,:,:) = Q(2:Nx,:,:)
        
        QR2(Nx - 1:Nx,:,:) = Q(Nx - 1:Nx,:,:)
        QR2(1:Nx - 2,:,:) = Q(3:Nx,:,:)
        
        QR3(Nx - 2:Nx,:,:) = Q(Nx - 2:Nx,:,:)
        QR3(1:Nx - 3,:,:) = Q(4:Nx,:,:)
        
    end if
    
    do i = 1,Nx
        do j = 1,Ny
            ! Step 3£ºcalculate the max eigenvalue
            rho = Q(i,j,1)
            u = Q(i,j,2)/Q(i,j,1)
            v = Q(i,j,3)/Q(i,j,1)
            E = Q(i,j,4)
            B1 = Q(i,j,5)
            B2 = Q(i,j,6)
    
            call wavespeeds(rho,u,v,E,B1,B2,1d0,0d0,AlphaV)
    
            AAlpha(:,1) = AlphaV(:,1)
            AAlpha(:,2) = AlphaV(:,1)
            AAlpha(:,3) = AlphaV(:,1)
            AAlpha(:,4) = AlphaV(:,1)
            AAlpha(:,5) = AlphaV(:,1)
            
            ! Step 4£ºcalculate f+ and f-
            
            ! f+
            fQP(1:dim,1) = fQL2(i,j,:)
            fQP(1:dim,2) = fQL1(i,j,:)
            fQP(1:dim,3) = fQ(i,j,:)
            fQP(1:dim,4) = fQR1(i,j,:)
            fQP(1:dim,5) = fQR2(i,j,:)
            
            QP(1:dim,1) = QL2(i,j,:)
            QP(1:dim,2) = QL1(i,j,:)
            QP(1:dim,3) = Q(i,j,:)
            QP(1:dim,4) = QR1(i,j,:)
            QP(1:dim,5) = QR2(i,j,:)
            
            ! f-
            fQN(1:dim,1) = fQR3(i,j,:)
            fQN(1:dim,2) = fQR2(i,j,:)
            fQN(1:dim,3) = fQR1(i,j,:)
            fQN(1:dim,4) = fQ(i,j,:)
            fQN(1:dim,5) = fQL1(i,j,:)
            
            QN(1:dim,1) = QR3(i,j,:)
            QN(1:dim,2) = QR2(i,j,:)
            QN(1:dim,3) = QR1(i,j,:)
            QN(1:dim,4) = Q(i,j,:)
            QN(1:dim,5) = QL1(i,j,:)
            
            call dimtodim1(fQP)
            call dimtodim1(fQN)
            call dimtodim1(QP)
            call dimtodim1(QN)
            
            ! Step 5£ºcalculate the eigenmatrix
            rho = 0.5d0*(Q(i,j,1) + QR1(i,j,1))
            u = 0.5d0*(Q(i,j,2) + QR1(i,j,2))/rho
            v = 0.5d0*(Q(i,j,3) + QR1(i,j,3))/rho
            E = 0.5d0*(Q(i,j,4) + QR1(i,j,4))
            B1 = 0.5d0*(Q(i,j,5) + QR1(i,j,5))
            B2 = 0.5d0*(Q(i,j,6) + QR1(i,j,6))
            
            call eigenmatrix(R,L,rho,u,v,0d0,E,B1,B2,0d0,1d0,0d0,0d0)
            
            if (undoeig == 1) then
                L = 0
                R = 0
                do d = 1,dim1
                    L(d,d) = 1
                    R(d,d) = 1
                end do
            end if
            
            fQP = matmul(L,fQP)
            fQN = matmul(L,fQN)
            QP = matmul(L,QP)
            QN = matmul(L,QN)
            fQP = 0.5d0*(fQP + AAlpha*QP)
            fQN = 0.5d0*(fQN - AAlpha*QN)
            
            ! Step 6£ºWENO reconstruction
            call WENO5V(QRnew(:,1),fQP(:,1),fQP(:,2),fQP(:,3),fQP(:,4),fQP(:,5))
            call WENO5V(QLnew(:,1),fQN(:,1),fQN(:,2),fQN(:,3),fQN(:,4),fQN(:,5))
            
            
            ! Step 7£ºproject back and get the reconstruct value at i+1/2,j
            QRnew = matmul(R,QRnew + QLnew)
            
            call dim1todim(QRnew)
            
            fhatx(i + 1,j,:) = QRnew(1:dim,1)
        end do
    end do
    
    if (boundcondition == 1) then
        fhatx(1,:,:) = fhatx(Nx + 1,:,:)
    elseif (boundcondition == 2) then
        fhatx(1,:,:) = fhatx(2,:,:)
    end if
    
    
    ! y-direction
    boundcondition = 1
    
    ! Step 1
    do i = 1,Nx
        do j = 1,Ny
            call f2(Q(i,j,:),yv)
            fQ(i,j,:) = yv
        end do
    end do
    
    ! Step 2
    
    if (boundcondition == 1) then 
        fQL2(:,1:2,:) = fQ(:,Ny - 1:Ny,:)
        fQL2(:,3:Ny,:) = fQ(:,1:Ny - 2,:)
        
        fQL1(:,1,:) = fQ(:,Ny,:)
        fQL1(:,2:Ny,:) = fQ(:,1:Ny - 1,:)
        
        fQR1(:,Ny,:) = fQ(:,1,:)
        fQR1(:,1:Ny - 1,:) = fQ(:,2:Ny,:)
        
        fQR2(:,Ny - 1:Ny,:) = fQ(:,1:2,:)
        fQR2(:,1:Ny - 2,:) = fQ(:,3:Ny,:)
        
        fQR3(:,Ny - 2:Ny,:) = fQ(:,1:3,:)
        fQR3(:,1:Ny - 3,:) = fQ(:,4:Ny,:)
        
        QL2(:,1:2,:) = Q(:,Ny - 1:Ny,:)
        QL2(:,3:Ny,:) = Q(:,1:Ny - 2,:)
        
        QL1(:,1,:) = Q(:,Ny,:)
        QL1(:,2:Ny,:) = Q(:,1:Ny - 1,:)
        
        QR1(:,Ny,:) = Q(:,1,:)
        QR1(:,1:Ny - 1,:) = Q(:,2:Ny,:)
        
        QR2(:,Ny - 1:Ny,:) = Q(:,1:2,:)
        QR2(:,1:Ny - 2,:) = Q(:,3:Ny,:)
        
        QR3(:,Ny - 2:Ny,:) = Q(:,1:3,:)
        QR3(:,1:Ny - 3,:) = Q(:,4:Ny,:)
        
    else if (boundcondition == 2) then
        
        fQL2(:,1:2,:) = fQ(:,1:2,:)
        fQL2(:,3:Ny,:) = fQ(:,1:Ny - 2,:)
        
        fQL1(:,1,:) = fQ(:,1,:)
        fQL1(:,2:Ny,:) = fQ(:,1:Ny - 1,:)
        
        fQR1(:,Ny,:) = fQ(:,Ny,:)
        fQR1(:,1:Ny - 1,:) = fQ(:,2:Ny,:)
        
        fQR2(:,Ny - 1:Ny,:) = fQ(:,Ny - 1:Ny,:)
        fQR2(:,1:Ny - 2,:) = fQ(:,3:Ny,:)
        
        fQR3(:,Ny - 2:Ny,:) = fQ(:,Ny - 2:Ny,:)
        fQR3(:,1:Ny - 3,:) = fQ(:,4:Ny,:)
        
        QL2(:,1:2,:) = Q(:,1:2,:)
        QL2(:,3:Ny,:) = Q(:,1:Ny - 2,:)
        
        QL1(:,1,:) = Q(:,1,:)
        QL1(:,2:Ny,:) = Q(:,1:Ny - 1,:)
        
        QR1(:,Ny,:) = Q(:,Ny,:)
        QR1(:,1:Ny - 1,:) = Q(:,2:Ny,:)
        
        QR2(:,Ny - 1:Ny,:) = Q(:,Ny - 1:Ny,:)
        QR2(:,1:Ny - 2,:) = Q(:,3:Ny,:)
        
        QR3(:,Ny - 2:Ny,:) = Q(:,Ny - 2:Ny,:)
        QR3(:,1:Ny - 3,:) = Q(:,4:Ny,:)
        
    end if
    
    
    do i = 1,Nx
        do j = 1,Ny
            
            ! Step 3
            rho = Q(i,j,1)
            u = Q(i,j,2)/Q(i,j,1)
            v = Q(i,j,3)/Q(i,j,1)
            E = Q(i,j,4)
            B1 = Q(i,j,5)
            B2 = Q(i,j,6)
    
            call wavespeeds(rho,u,v,E,B1,B2,0d0,1d0,AlphaV)
    
            AAlpha(:,1) = AlphaV(:,1)
            AAlpha(:,2) = AlphaV(:,1)
            AAlpha(:,3) = AlphaV(:,1)
            AAlpha(:,4) = AlphaV(:,1)
            AAlpha(:,5) = AlphaV(:,1)
            
            ! Step 4
            
            ! f+
            fQP(1:dim,1) = fQL2(i,j,:)
            fQP(1:dim,2) = fQL1(i,j,:)
            fQP(1:dim,3) = fQ(i,j,:)
            fQP(1:dim,4) = fQR1(i,j,:)
            fQP(1:dim,5) = fQR2(i,j,:)
            
            QP(1:dim,1) = QL2(i,j,:)
            QP(1:dim,2) = QL1(i,j,:)
            QP(1:dim,3) = Q(i,j,:)
            QP(1:dim,4) = QR1(i,j,:)
            QP(1:dim,5) = QR2(i,j,:)
            
            ! f-
            fQN(1:dim,1) = fQR3(i,j,:)
            fQN(1:dim,2) = fQR2(i,j,:)
            fQN(1:dim,3) = fQR1(i,j,:)
            fQN(1:dim,4) = fQ(i,j,:)
            fQN(1:dim,5) = fQL1(i,j,:)
            
            QN(1:dim,1) = QR3(i,j,:)
            QN(1:dim,2) = QR2(i,j,:)
            QN(1:dim,3) = QR1(i,j,:)
            QN(1:dim,4) = Q(i,j,:)
            QN(1:dim,5) = QL1(i,j,:)
            
            call dimtodim1(fQP)
            call dimtodim1(fQN)
            call dimtodim1(QP)
            call dimtodim1(QN)
            
            ! Step 5
            rho = 0.5d0*(Q(i,j,1) + QR1(i,j,1))
            u = 0.5d0*(Q(i,j,2) + QR1(i,j,2))/rho
            v = 0.5d0*(Q(i,j,3) + QR1(i,j,3))/rho
            E = 0.5d0*(Q(i,j,4) + QR1(i,j,4))
            B1 = 0.5d0*(Q(i,j,5) + QR1(i,j,5))
            B2 = 0.5d0*(Q(i,j,6) + QR1(i,j,6))
            
            call eigenmatrix(R,L,rho,u,v,0d0,E,B1,B2,0d0,0d0,1d0,0d0)
            
            if (undoeig == 1) then
                L = 0
                R = 0
                do d = 1,dim1
                    L(d,d) = 1
                    R(d,d) = 1
                end do
            end if
            
            fQP = matmul(L,fQP)
            fQN = matmul(L,fQN)
            QP = matmul(L,QP)
            QN = matmul(L,QN)
            fQP = 0.5d0*(fQP + AAlpha*QP)
            fQN = 0.5d0*(fQN - AAlpha*QN)
            
            ! Step 6
            call WENO5V(QRnew(:,1),fQP(:,1),fQP(:,2),fQP(:,3),fQP(:,4),fQP(:,5))
            call WENO5V(QLnew(:,1),fQN(:,1),fQN(:,2),fQN(:,3),fQN(:,4),fQN(:,5))
            
            
            ! Step 7
            QRnew = matmul(R,QRnew + QLnew)
            
            call dim1todim(QRnew)
            
            fhaty(i,j + 1,:) = QRnew(1:dim,1)
            
        end do
    end do
    
    if (boundcondition == 1) then
        fhaty(:,1,:) = fhaty(:,Ny + 1,:)
    elseif (boundcondition == 2) then
        fhaty(:,1,:) = fhaty(:,2,:)
    end if
    
    DQ = -(fhatx(2:Nx + 1,:,:) - fhatx(1:Nx,:,:))/hx - (fhaty(:,2:Ny + 1,:) - fhaty(:,1:Ny,:))/hy
    
    end subroutine Lh
    
    
    
    
    
    
    
    
    
    subroutine WENO5(h,a,b,c,d,e)
    
    real h,a,b,c,d,e,epsilon,h0,h1,h2,beta0,beta1,beta2,w0,w1,w2,S
    
    epsilon = 1e-6;
    
    h0 = (2*a - 7*b + 11*c)/6d0
    h1 = (-b + 5*c + 2*d)/6d0
    h2 = (2*c + 5*d - e)/6d0
    
    beta0 = (13d0/12d0)*(a - 2*b + c)**2 + 0.25d0*(a - 4*b + 3*c)**2
    beta1 = (13d0/12d0)*(b - 2*c + d)**2 + 0.25d0*(b - d)**2
    beta2 = (13d0/12d0)*(c - 2*d + e)**2 + 0.25d0*(3*c - 4*d + e)**2
    
    w0 = 1d0/(epsilon + beta0)**2
    w1 = 6d0/(epsilon + beta1)**2
    w2 = 3d0/(epsilon + beta2)**2
    S = w0 + w1 + w2
    w0 = w0/S
    w1 = w1/S
    w2 = w2/S
    
    h = w0*h0 + w1*h1 + w2*h2
    
    end subroutine WENO5
    
    
    
    
    subroutine WENO5V(y,a,b,c,d,e)
    
    integer dim
    parameter(dim = 8)
    real y(dim),a(dim),b(dim),c(dim),d(dim),e(dim),yk
    integer k
    
    do k = 1,dim
        call WENO5(yk,a(k),b(k),c(k),d(k),e(k))
        y(k) = yk
    end do
    
    end subroutine WENO5V
    
    subroutine dimtodim1(x)
    
    real x(8,5)
    
    x(6:7,:) = x(5:6,:)
    x(5,:) = x(4,:)
    x(4,:) = 0
    x(8,:) = 0
    
    end subroutine dimtodim1
    
    subroutine dim1todim(x)
    
    real x(8,1)
    
    x(4,:) = x(5,:)
    x(5:6,:) = x(6:7,:)
    x(7,:) = 0
    x(8,:) = 0
    
    end subroutine dim1todim
    
    
    subroutine wavespeeds(rho,u1,u2,E,B1,B2,n1,n2,alpha)
    ! max wave speed, 2D Euler equation.  n1**2 + n2**2 need not == 1
    implicit none
    real :: rho, u1, u2, E, B1, B2, n1, n2, gamma,c,ca,cf,cs,un,p
    real :: alpha(8,1)
    
    gamma = 5d0/3d0
    
    un = u1*n1 + u2*n2
    
    p = (gamma - 1)*(E - 0.5d0*rho*(u1**2 + u2**2) - 0.5d0*(B1**2 + B2**2))
    
    c = (abs(gamma*p/rho))**0.5d0
    
    ca = (abs((B1*n1 + B2*n2)**2/rho))**0.5d0
    
    cf = abs( 0.5d0*( c**2 + (B1**2 + B2**2)/rho + ( (abs(c**2 + (B1**2 + B2**2)/rho)**2 -4*c**2*(B1*n1 + B2*n2)**2/rho ))**0.5d0 ) )**0.5d0
    
    cs = abs( 0.5d0*( c**2 + (B1**2 + B2**2)/rho - ( (abs(c**2 + (B1**2 + B2**2)/rho)**2 -4*c**2*(B1*n1 + B2*n2)**2/rho ))**0.5d0 ) )**0.5d0
    
    alpha(1,1) = un
    alpha(2,1) = un
    alpha(3,1) = un + ca
    alpha(4,1) = un - ca
    alpha(5,1) = un + cf
    alpha(6,1) = un - cf
    alpha(7,1) = un + cs
    alpha(8,1) = un - cs
    
    alpha = abs(alpha)
    
  end subroutine wavespeeds
    
    