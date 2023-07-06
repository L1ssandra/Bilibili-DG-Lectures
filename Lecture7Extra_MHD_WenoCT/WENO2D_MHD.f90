    ! WENO2D_MHD.f90 
    
    ! The WENO-CT method for 2D ideal MHD equation
    program main

    implicit none

    ! variables

    real pi,CFL,x,y,t,xa,xb,ya,yb,dt,tend,gamma,hx,hy
    integer i,j,k,n,Nx,Ny,RKorder,dim,d,frame,skiptime,Nx1,Ny1
    integer Qlength,countstep,count,skip
    real stoptime1,stoptime2,stoptime3,umax
    
    ! The parameter of each example shoule be:
    
    ! Omega = [-10,10]*[-10,10], Tend = 2 for SMOOTH VORTEX
    ! Omega = [0,2*pi]*[0,2*pi], Tend = 4 for ORSZAG-TANG VORTEX
    ! Omega = [0,1]*[0,1], Tend = 0.295 for ROTOR

    parameter(dim = 6)
    parameter(pi = 4*atan(1.0d0))
    parameter(CFL = 3d0)
    parameter(Nx = 400)
    parameter(Ny = 400)
    parameter(Nx1 = Nx + 1)
    parameter(Ny1 = Ny + 1)
    parameter(xa = 0)
    parameter(xb = 1)
    parameter(ya = 0)
    parameter(yb = 1)
    parameter(hx = (xb - xa)/Nx)
    parameter(hy = (yb - ya)/Ny)
    parameter(tend = 0.295)
    parameter(RKorder = 4)
    parameter(gamma = 5d0/3d0)
    parameter(frame = 5) ! the max frame
    parameter(Qlength = Nx*Ny*frame)
    parameter(skiptime = 5000) ! save the solution at t = skiptime*N
    parameter(stoptime1 = 0.5d0)
    parameter(stoptime2 = 2.0d0)
    parameter(stoptime3 = 3.0d0)
    
    real XX(Nx + 1),YY(Ny + 1),Xc(Nx),Yc(Ny),Q(Nx,Ny,dim),QR(Nx,Ny,dim),Linfty(dim),Linfty1,L2,L2V(dim),A(Nx,Ny),QA(Nx,Ny,dim),B(Nx,Ny,2)
    real TT(frame),QF(Nx,Ny,dim,frame),Error(Nx,Ny,dim),AR(Nx,Ny)
    real QQ1(1,Qlength),QQ2(1,Qlength),QQ3(1,Qlength),QQ4(1,Qlength),QQ5(1,Qlength),QQ6(1,Qlength),AF(Nx,Ny,frame),AA(1,Qlength)
    real alphax,alphay,alpha1,alpha2
    integer :: countstop = 0
    real,external:: eigxmax,eigymax,rho0,ux0,uy0
    
    ! The boundary conditions (boundaryX,boundaryY) should be:
    
    ! SMOOTH VORTEX: (1,1)
    ! ORSZAG-TANG VORTEX: (1,1).
    ! ROTOR: (1,1) in "RK4Chara" and (1,3) in "FDcurl"&"RK4HJ".
    
    ! initial data
    real(kind = 8) r2
    r2(x,y) = x**2.0d0 + y**2.0d0
    real p
    !p(x,y) = 1 - r2(x,y)/(8.0d0*pi**2)*exp(1 - r2(x,y)) ! smooth vortex
    !p(x,y) = gamma ! Orszag-Tang Vortex
    p(x,y) = 0.5d0 ! Rotor
    real rho
    !rho(x,y) = 1 ! smooth vortex
    !rho(x,y) = gamma**2 ! Orszag-Tang Vortex
    rho(x,y) = rho0(x,y) ! Rotor
    real v1
    !v1(x,y) = 1 - 1d0/(2.0d0*pi)*y*exp(0.5d0*(1 - r2(x,y))) ! smooth vortex
    !v1(x,y) = -sin(y) ! Orszag-Tang Vortex
    v1(x,y) = ux0(x,y) ! Rotor
    real v2
    !v2(x,y) = 1 + 1d0/(2.0d0*pi)*x*exp(0.5d0*(1 - r2(x,y))) ! smooth vortex
    !v2(x,y) = sin(x) ! Orszag-Tang Vortex
    v2(x,y) = uy0(x,y) ! Rotor
    real B1
    !B1(x,y) = -1d0/(2.0d0*pi)*y*exp(0.5d0*(1 - r2(x,y))) ! smooth vortex
    !B1(x,y) = -sin(y) ! Orszag-Tang Vortex
    B1(x,y) = 2.5d0/(4d0*pi)**0.5d0 ! Rotor
    real B2
    !B2(x,y) = 1d0/(2.0d0*pi)*x*exp(0.5d0*(1 - r2(x,y))) ! smooth vortex
    !B2(x,y) = sin(2d0*x) ! Orszag-Tang Vortex
    B2(x,y) = 0 ! Rotor
    real U1
    U1(x,y) = rho(x,y)
    real U2
    U2(x,y) = rho(x,y)*v1(x,y)
    real U3
    U3(x,y) = rho(x,y)*v2(x,y)
    real U4
    U4(x,y) = p(x,y)/(gamma - 1) + 0.5d0*rho(x,y)*(v1(x,y)**2 + v2(x,y)**2) + 0.5d0*(B1(x,y)**2 + B2(x,y)**2)
    real U5
    U5(x,y) = B1(x,y)
    real U6
    U6(x,y) = B2(x,y)
    
    real A0
    !A0(x,y) = 1d0/(2d0*pi)*exp(0.5d0*(1 - r2(x,y))) ! smooth vortex
    !A0(x,y) = 0.5d0*cos(2d0*x) + cos(y) ! Orszag-Tang Vortex
    A0(x,y) = 2.5d0*y/(4d0*pi)**0.5d0 ! Rotor
    

    do i = 1,Nx + 1
        XX(i) = xa + hx*(i - 1)
    end do

    do j = 1,Ny + 1
        YY(j) = ya + hy*(j - 1)
    end do

    Xc = 0.5d0*(XX(1:Nx) + XX(2:Nx + 1))
    Yc = 0.5d0*(YY(1:Ny) + YY(2:Ny + 1))

    ! The value at t = 0
    do i = 1,Nx
        do j = 1,Ny
            Q(i,j,1) = U1(Xc(i),Yc(j))
            Q(i,j,2) = U2(Xc(i),Yc(j))
            Q(i,j,3) = U3(Xc(i),Yc(j))
            Q(i,j,4) = U4(Xc(i),Yc(j))
            Q(i,j,5) = U5(Xc(i),Yc(j))
            Q(i,j,6) = U6(Xc(i),Yc(j))
            A(i,j) = A0(Xc(i),Yc(j))
        end do
    end do

    ! real solution
    QR = Q
    AR = A
    ! flash
    QF(:,:,:,1) = Q
    AF(:,:,1) = A

    t = 0
    TT(1) = t
    count = 2
    skip = 0

    ! start
    do while (t < tend)

        alphax = 1
        alphay = 1

        do i = 1,Nx
            do j = 1,Ny
                alpha1 = eigxmax(Q(i,j,:))
                alpha2 = eigymax(Q(i,j,:))
                if (alpha1 > alphax) then
                    alphax = alpha1
                end if
                if (alpha2 > alphay) then
                    alphay = alpha2
                end if
            end do
        end do

        dt = CFL/(alphax/hx + alphay/hy)

        if (t + dt <= tend) then
            t = t + dt
        else
            dt = tend - t
            t = tend
        end if
        
        if (dt <= 1e-6) then
            dt = tend - t
            t = tend
        end if
        
        ! save the solution at some times
        if ((t >= stoptime1) .and. (countstop == 0)) then
            dt = stoptime1 + dt - t
            t = stoptime1
            countstop = countstop + 1
            skip = skiptime - 1
        end if
        
        if ((t >= stoptime2) .and. (countstop == 1)) then
            dt = stoptime2 + dt - t
            t = stoptime2
            countstop = countstop + 1
            skip = skiptime - 1
        end if
        
        if ((t >= stoptime3) .and. (countstop == 2)) then
            dt = stoptime3 + dt - t
            t = stoptime3
            countstop = countstop + 1
            skip = skiptime - 1
        end if
        
        if (t == tend) then
            skip = skiptime - 1
        end if
        
        skip = skip + 1
        
        
        
        if (skip == skiptime) then
            skip = 0
        end if
        
        umax = 0
        ! calculate_umax
        do i = 1,Nx
            do j = 1,Ny
                
                if (Q(i,j,1) > umax) then
                    umax = Q(i,j,1)
                end if
                
            end do
        end do
        
        print *,"t = ",t,"umax = ",umax
        
        if (RKorder == 4) then
            call RK4Chara(Q,A,hx,hy,dt,Nx,Ny)
        end if
        
        if (skip == 0) then
    
            TT(count) = t
            QF(:,:,:,count) = Q
            AF(:,:,count) = A
            print *,"save the solution at t = ",t
            count = count + 1
            
        end if
     
    end do

    ! calculate error
    Error = abs(Q - QR)
    Linfty = 0
    L2V = 0

    do d = 1,dim
        L2 = 0
        do i = 1,Nx
            do j = 1,Ny
                L2 = L2 + Error(i,j,d)**2
                if (Error(i,j,d) > Linfty(d)) then
                    Linfty(d) = Error(i,j,d)
                end if
            end do
        end do
        L2V(d) = (L2/(Nx*Ny))**0.5d0
    end do
    
    print *,L2V
    print *,Linfty

    QQ1 = reshape(QF(:,:,1,:),(/1,Qlength/))
    QQ2 = reshape(QF(:,:,2,:),(/1,Qlength/))
    QQ3 = reshape(QF(:,:,3,:),(/1,Qlength/))
    QQ4 = reshape(QF(:,:,4,:),(/1,Qlength/))
    QQ5 = reshape(QF(:,:,5,:),(/1,Qlength/))
    QQ6 = reshape(QF(:,:,6,:),(/1,Qlength/))
    AA = reshape(AF,(/1,Qlength/))
    
    open(unit = 2,file = 'T.txt')
        do i = 1,count - 1
            write(2,*) TT(i)
        end do
        
    print *,"writing"
    
    open(unit = 1,file = 'Q1.txt')
    open(unit = 3,file = 'Q2.txt')
    open(unit = 4,file = 'Q3.txt')
    open(unit = 5,file = 'Q4.txt')
    open(unit = 6,file = 'Q5.txt')
    open(unit = 7,file = 'Q6.txt')
    open(unit = 8,file = 'A.txt')
        do i = 1,Nx*Ny*(count - 1)
            if (i <= Nx*Ny*count) then
                write(1,*) QQ1(1,i)
                write(3,*) QQ2(1,i)
                write(4,*) QQ3(1,i)
                write(5,*) QQ4(1,i)
                write(6,*) QQ5(1,i)
                write(7,*) QQ6(1,i)
                write(8,*) AA(1,i)
                if (mod(i,5000) == 0) then
                    print *,i,"/",Nx*Ny*(count - 1)
                end if
            end if
        end do
        
    open(unit = 9,file = 'Nx.txt')
        write(9,*) Nx
    open(unit = 10,file = 'Ny.txt')
        write(10,*) Ny
    
    
    
    

    end program main

    subroutine f(x,y)


    real x(6),y(6,2),gamma,v1,v2,p,rho,S,T,K

    gamma = 5d0/3d0

    p = (x(4) - 0.5d0*(x(2)**2 + x(3)**2)/x(1) - 0.5d0*(x(5)**2 + x(6)**2))*(gamma - 1)
    S = p + 0.5d0*(x(5)**2 + x(6)**2)
    T = x(4) + p + 0.5d0*(x(5)**2 + x(6)**2)
    K = x(2)*x(5)/x(1) + x(3)*x(6)/x(1)
    
    y(1,1) = x(2)                                                    
    y(1,2) = x(3)
    y(2,1) = x(2)**2/x(1) + S - x(5)**2
    y(2,2) = x(2)*x(3)/x(1) - x(5)*x(6)
    y(3,1) = x(2)*x(3)/x(1) - x(5)*x(6)                            
    y(3,2) = x(3)**2/x(1) + S - x(6)**2
    y(4,1) = T*x(2)/x(1) - K*x(5)
    y(4,2) = T*x(3)/x(1) - K*x(6)
    y(5,1) = 0
    y(5,2) = -(x(2)*x(6) - x(3)*x(5))/x(1)
    y(6,1) = -(x(3)*x(5) - x(2)*x(6))/x(1)
    y(6,2) = 0

    end subroutine f
    
    
    
    
    
    subroutine f1(u,y)
    
    real u(6),y(6),y2(6,2)
    
    call f(u,y2)
    
    y = y2(:,1)
    
    end subroutine f1
    
    
    
    
    
    subroutine f2(u,y)
    
    real u(6),y(6),y2(6,2)
    
    call f(u,y2)
    
    y = y2(:,2)
    
    end subroutine f2
    
    
    ! 求最大特征值
    function eigxmax(xv)
    
    real xv(6),a,cf
    real :: eigxmax
    
    a = abs((5.0d0/3.0d0)*(xv(4) - 0.5d0*(xv(2)**2.0d0 + xv(3)**2.0d0)/xv(1) - 0.5d0*(xv(5)**2 + xv(6)**2))/xv(1))**0.5d0
    cf = 0.5d0*abs(a**2 + (xv(5)**2 + xv(6)**2)/xv(1) + ((a**2 + (xv(5)**2 + xv(6)**2)/xv(1))**2 - 4*a**2*(xv(5))**2/xv(1))**0.5d0 )**0.5d0
    eigxmax = abs(xv(2)/xv(1)) + abs(cf)
    !eigxmax = 1
    
    return
    
    end
    
    
    
    
    function eigymax(xv)
    
    implicit none
    
    real xv(6),a,cf
    real :: eigymax
    
    a = abs((5.0d0/3.0d0)*(xv(4) - 0.5d0*(xv(2)**2.0d0 + xv(3)**2.0d0)/xv(1) - 0.5d0*(xv(5)**2 + xv(6)**2))/xv(1))**0.5d0
    cf = 0.5d0*abs(a**2 + (xv(5)**2 + xv(6)**2)/xv(1) + ((a**2 + (xv(5)**2 + xv(6)**2)/xv(1))**2 - 4*a**2*(xv(6))**2/xv(1))**0.5d0 )**0.5d0
    eigymax = abs(xv(3)/xv(1)) + abs(cf)
    !eigymax = 1
    
    return
    
    end
    
    
    function rho0(x,y)
    
    implicit none
    
    real x,y,f,r,r0,r1
    real rho0
    
    r = ((x - 0.5d0)**2d0 + (y - 0.5d0)**2d0)**0.5d0
    r0 = 0.1d0
    r1 = 0.115d0
    
    f = (r1 - r)/(r1 - r0)
    
    if (r < r0) then
        rho0 = 10d0
    else if ((r >= r0) .and. (r < r1)) then
        rho0 = 1d0 + 9d0*f
    else if (r >= r1) then
        rho0 = 1d0
    end if
    
    end
    
    
    
    function ux0(x,y)
    
    implicit none
    
    real x,y,f,r,r0,r1
    real ux0
    
    r = ((x - 0.5d0)**2d0 + (y - 0.5d0)**2d0)**0.5d0
    r0 = 0.1d0
    r1 = 0.115d0
    
    f = (r1 - r)/(r1 - r0)
    
    if (r < r0) then
        ux0 = -(y - 0.5d0)/r0
    else if ((r >= r0) .and. (r < r1)) then
        ux0 = -f*(y - 0.5d0)/r
    else if (r >= r1) then
        ux0 = 0d0        
    end if
    
    end
    
    
    
    function uy0(x,y)
    
    implicit none
    
    real x,y,f,r,r0,r1
    real uy0
    
    r = ((x - 0.5d0)**2d0 + (y - 0.5d0)**2d0)**0.5d0
    r0 = 0.1d0
    r1 = 0.115d0
    
    f = (r1 - r)/(r1 - r0)
    
    if (r < r0) then
        uy0 = (x - 0.5d0)/r0
    else if ((r >= r0) .and. (r < r1)) then
        uy0 = f*(x - 0.5d0)/r
    else if (r >= r1) then
        uy0 = 0d0        
    end if
    
    end
    
    
    

