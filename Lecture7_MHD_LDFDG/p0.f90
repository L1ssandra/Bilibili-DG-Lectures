    function p0(x,y)
    
    real x,y,r,r0,p0
    
    r = (x**2 + y**2)**0.5
    r0 = 0.1
    
    if (r <= r0) then
        p0 = 1000
    else
        p0 = 0.1
    end if
    
    end