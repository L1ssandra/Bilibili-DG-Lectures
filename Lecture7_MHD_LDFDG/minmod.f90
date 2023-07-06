    subroutine minmod
    
    include 'com.txt'
    
    if (direction == 1) then
        hd = hx
    else if (direction == 2) then
        hd = hy
    end if
    
    do i = 1,8
        if (abs(DeltaU(i,1)) <= M*hd**2) then
            DeltaUmod(i,1) = DeltaU(i,1)
        else
            a = sign(1d0,DeltaU(i,1))
            b = sign(1d0,DeltaUR(i,1))
            c = sign(1d0,DeltaUL(i,1))
            s = (a + b + c)/3d0
            if (abs(s) == 1) then
                !DeltaUmod(i,1) = s*min(abs(DeltaU(i,1)),0.5*abs(DeltaUR(i,1)),0.5*abs(DeltaUL(i,1)))
                DeltaUmod(i,1) = s*min(abs(DeltaU(i,1)),abs(DeltaUR(i,1)),abs(DeltaUL(i,1)))
            else
                DeltaUmod(i,1) = 0
            end if
        end if
        
    end do
    
    end subroutine minmod