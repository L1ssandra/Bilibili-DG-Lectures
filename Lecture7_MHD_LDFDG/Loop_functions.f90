    function Bx01(x,y)
    
    real x,y,A0,r0
    real Bx01
    
    r = (x**2 + y**2)**0.5
    r0 = 0.3
    A0 = 1e-3
    
    if (r < r0) then
        Bx01 = -A0*y/r
    else
        Bx01 = 0
    end if
    
    end
    
    
    function By01(x,y)
    
    real x,y,A0,r0
    real By01
    
    r = (x**2 + y**2)**0.5
    r0 = 0.3
    A0 = 1e-3
    
    if (r < r0) then
        By01 = A0*x/r
    else
        By01 = 0
    end if
    
    end
    
    function Az(x,y)
    
    real x,y,A0,r0
    real Az
    
    r = (x**2 + y**2)**0.5
    r0 = 0.3
    A0 = 1e-3
    
    if (r < r0) then
        Az = A0*(r0 - r)
    else
        Az = 0
    end if
    
    end