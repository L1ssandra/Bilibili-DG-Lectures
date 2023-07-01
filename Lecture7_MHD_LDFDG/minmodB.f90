    subroutine minmodB(a1,a,b,c,M,hd)
    
    real a1,a,b,c,as,bs,cs,M,hd
    
    if (abs(a) <= M*hd**2) then
        a1 = a
    else
        as = sign(1d0,a)
        bs = sign(1d0,b)
        cs = sign(1d0,c)
        s = (as + bs + cs)/3d0
        if (abs(s) == 1) then
            !a1 = s*min(abs(a),0.5*abs(b),0.5*abs(c))
            a1 = s*min(abs(a),abs(b),abs(c))
            !if (a1 == b .or. a1 == c) then
            !    a1 = 0.5*a1
            !end if
        else
            a1 = 0
        end if
    end if
    
    end subroutine minmodB