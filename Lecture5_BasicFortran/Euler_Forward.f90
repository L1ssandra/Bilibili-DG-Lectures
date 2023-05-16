    subroutine Euler_Forward
    
    include 'com.txt'
    
    CFL = 0.2
    dt = CFL*hx
    t = 0
    t1 = tend/frameMAX
    i1 = 1
    
    open(unit = 1,file = 'uh.txt')
    open(unit = 2,file = 'T.txt')
    
    do j = 0,Ny
        do i = 0,Nx
            write(1,*) uh(i,j)
        end do
    end do
    write(2,*) t
    
    do while (t < tend)
        
        if (t + dt > tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        call Lh
        
        uh = uh + dt*du
        
        call calculate_maxu
        
        print *,t,maxu
        
        if (t >= i1*t1) then
            do j = 0,Ny
                do i = 0,Nx
                    write(1,*) uh(i,j)
                end do
            end do
            write(2,*) t
            print *,"save the solution at t = ",t,i1*t1
            i1 = i1 + 1
        end if
        
    end do
    
    end subroutine Euler_Forward