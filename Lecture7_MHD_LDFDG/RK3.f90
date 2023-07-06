    subroutine RK3
    
    include 'com.txt'
    
    t = 0
    count = 0
    
    call set_bc
    
    call TVB_Limiter_P2
    
    !call div_free
    
    call calculate_umax
        
    call calculate_totaldiv
    
    call calculate_pmin
    
    !if (totaldiv > 1e-7) then
    !    call calculate_Az
    !    call div_free_Balsara
    !    call calculate_totaldiv
    !    call calculate_umax
    !end if
        
    print *,t,"  ",umax,"  ",totaldiv
    
    open(unit = 1,file = 'Q1flash.txt')
    open(unit = 2,file = 'Q2flash.txt')
    open(unit = 3,file = 'Q3flash.txt')
    open(unit = 4,file = 'Q4flash.txt')
    open(unit = 5,file = 'Q5flash.txt')
    open(unit = 6,file = 'Q6flash.txt')
    open(unit = 7,file = 'Q7flash.txt')
    open(unit = 8,file = 'Q8flash.txt')
    open(unit = 9,file = 'T.txt')
    
    if (flash == 1) then
        do d = 1,dimPk
            do j = 1,Ny
                do i = 1,Nx
                    write(1,*) uh(i,j,d,1)
                    write(2,*) uh(i,j,d,2)
                    write(3,*) uh(i,j,d,3)
                    write(4,*) uh(i,j,d,4)
                    write(5,*) uh(i,j,d,5)
                    write(6,*) uh(i,j,d,6)
                    write(7,*) uh(i,j,d,7)
                    write(8,*) uh(i,j,d,8)
                end do
            end do
        end do
        write(9,*) t
    end if
    
    tt = tend/frameMAX
    countf = 1
    
    do while (t < tend)
        
        call calculate_dt
    
        if (t + dt >= tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        if (dt < 1e-7) then
            t = tend
        end if
        
        ! Stage 1
        call set_bc
    
        call Lh
        
        call LhEz
        
        uI = uh + dt*du
        
        BxI = Bx + dt*dEz1
        ByI = By + dt*dEz2
        
        uh0 = uh
        Bx0 = Bx
        By0 = By
        
        uh = uI
        Bx = BxI
        By = ByI
        
        call TVB_Limiter_P2
        
        !call div_free
            
        call set_bc
        
        if (RKorder == 3) then
            
            !Stage 2
            call Lh
            
            call LhEz
            
            uh = uI + dt*du
            Bx = BxI + dt*dEz1
            By = ByI + dt*dEz2
        
            uII = (3d0/4d0)*uh0 + (1d0/4d0)*uh
            BxII = (3d0/4d0)*Bx0 + (1d0/4d0)*Bx
            ByII = (3d0/4d0)*By0 + (1d0/4d0)*By
        
            uh = uII
            Bx = BxII
            By = ByII
        
            call TVB_Limiter_P2
            
            !call div_free
            
            !Stage 3
            call set_bc
        
            call Lh
            
            call LhEz
            
            Bx = BxII + dt*dEz1
            By = ByII + dt*dEz2
            uh = uII + dt*du
        
            uh = (1d0/3d0)*uh0 + (2d0/3d0)*uh
            Bx = (1d0/3d0)*Bx0 + (2d0/3d0)*Bx
            By = (1d0/3d0)*By0 + (2d0/3d0)*By
        
            call TVB_Limiter_P2
            
            !call div_free
            
        end if
        
        count = count + 1
        
        call calculate_umax
        
        call calculate_totaldiv
        
        call calculate_pmin
        
        print *,t,"  ",umax,"  ",totaldiv
        
        ! save solution
        if ((t > tt*countf .or. t == tend) .and. flash == 1) then
            print *,"Ð´ÈëÊý¾Ý",t,tt*countf
            do d = 1,dimPk
                do j = 1,Ny
                    do i = 1,Nx
                        write(1,*) uh(i,j,d,1)
                        write(2,*) uh(i,j,d,2)
                        write(3,*) uh(i,j,d,3)
                        write(4,*) uh(i,j,d,4)
                        write(5,*) uh(i,j,d,5)
                        write(6,*) uh(i,j,d,6)
                        write(7,*) uh(i,j,d,7)
                        write(8,*) uh(i,j,d,8)
                    end do
                end do
            end do
            write(9,*) t
            countf = countf + 1
        end if
        
    end do
    
    end subroutine RK3