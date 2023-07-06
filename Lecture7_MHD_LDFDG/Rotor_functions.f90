    function rho0(x,y)
    
    implicit none
    
    real(kind = 8) x,y,f,r,r0,r1
    real(kind = 8) rho0
    
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
    
    real(kind = 8) x,y,f,r,r0,r1
    real(kind = 8) ux0
    
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
    
    real(kind = 8) x,y,f,r,r0,r1
    real(kind = 8) uy0
    
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