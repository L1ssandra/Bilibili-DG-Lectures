    function rhoMa(x,y)
    
    real x,y
    real rhoMa
    
    if (x < 1d0/6d0 + y/3d0**0.5) then
        rhoMa = 8
    else
        rhoMa = 1.4
    end if
    
    end
    !*********************************************
    function v1Ma(x,y)
    
    parameter(pi = 4*atan(1.0d0))
    real x,y
    real v1Ma
    
    if (x < 1d0/6d0 + y/3d0**0.5) then
        v1Ma = 8.25*cos(pi/6d0)
    else
        v1Ma = 0
    end if
    
    end
    !*********************************************
    function v2Ma(x,y)
    
    parameter(pi = 4*atan(1.0d0))
    real x,y
    real v2Ma
    
    if (x < 1d0/6d0 + y/3d0**0.5) then
        v2Ma = -8.25*sin(pi/6d0)
    else
        v2Ma = 0
    end if
    
    end
    !*********************************************
    function pMa(x,y)
    
    real x,y
    real pMa
    
    if (x < 1d0/6d0 + y/3d0**0.5) then
        pMa = 116.5
    else
        pMa = 1
    end if
    
    end