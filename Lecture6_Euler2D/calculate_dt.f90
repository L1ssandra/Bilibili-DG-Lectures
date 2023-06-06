    subroutine calculate_dt
    
    include 'com.txt'
    
    if (RKorder == 1) then
        CFL = 0.05
    else if (RKorder == 3) then
        CFL = 0.2
    end if
    
    alphax = 0
    alphay = 0
    
    do i = 0,Nx1
        do j = 0,Ny1
            call eigenvalueMm(alpha1,alpha2,uh(i,j,1,1),uh(i,j,1,2),uh(i,j,1,3),uh(i,j,1,4),1,0)
            if (abs(alpha1) > alphax .or. abs(alpha2) > alphax) then
                alphax = max(abs(alpha1),abs(alpha2))
            end if
            call eigenvalueMm(alpha1,alpha2,uh(i,j,1,1),uh(i,j,1,2),uh(i,j,1,3),uh(i,j,1,4),0,1)
            if (abs(alpha1) > alphay .or. abs(alpha2) > alphay) then
                alphay = max(abs(alpha1),abs(alpha2))
            end if
        end do
    end do
    
    dt = CFL/(alphax/hx + alphay/hy)
    
    end subroutine calculate_dt