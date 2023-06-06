    subroutine RK3
    
    include 'com.txt'
    
    t = 0
    tRK = 0
    count = 0
    
    call set_bc
    
    call TVB_Limiter_P2
    
    call calculate_umax
    
    call calculate_pmin
        
    print *,t,"  ",umax,"  ",pmin
    
    open(unit = 1,file = 'Q1flash.txt')
    open(unit = 2,file = 'Q2flash.txt')
    open(unit = 3,file = 'Q3flash.txt')
    open(unit = 4,file = 'Q4flash.txt')
    open(unit = 9,file = 'T.txt')
    
    do d = 1,dimPk
        do j = 1,Ny
            do i = 1,Nx
                write(1,*) uh(i,j,d,1)
                write(2,*) uh(i,j,d,2)
                write(3,*) uh(i,j,d,3)
                write(4,*) uh(i,j,d,4)
            end do
        end do
    end do
    write(9,*) t
    
    tt = tend/frameMAX
    countf = 1
    
    do while (t < tend)
        
        call calculate_dt
        
        tRK = t
    
        if (t + dt >= tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        if (dt < 1e-8) then
            t = tend
        end if
        
        ! Stage 1
        call set_bc
    
        call Lh
        
        uI = uh + dt*du
        
        uh0 = uh
        
        uh = uI
        
        call TVB_Limiter_P2
        
        call pp_Limiter
        
        if (RKorder == 3) then
            
            !Stage 2
            tRK = tRK + dt
            
            call set_bc
            
            call Lh
            
            uh = uI + dt*du
        
            uII = (3d0/4d0)*uh0 + (1d0/4d0)*uh
        
            uh = uII
        
            call TVB_Limiter_P2
            
            call pp_Limiter
            
            !Stage 3
            tRK = tRK - 0.5*dt
            
            call set_bc
        
            call Lh
            
            uh = uII + dt*du
        
            uh = (1d0/3d0)*uh0 + (2d0/3d0)*uh
        
            call TVB_Limiter_P2
            
            call pp_Limiter
            
        end if
        
        count = count + 1
        
        call calculate_umax
        
        call calculate_pmin
        
        print *,t,"  ",umax,"  ",pmin
        
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
                    end do
                end do
            end do
            write(9,*) t
            countf = countf + 1
        end if
        
    end do
    
    end subroutine RK3