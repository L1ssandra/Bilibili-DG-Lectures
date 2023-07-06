    function rhoC2(x,y)
    
    real x,y,r,r0
    real rhoC2
    
    r = ((x - 0.25)**2 + (y - 0.5)**2)**0.5
    r0 = 0.15
    
    if (x <= 0.05) then ! Omega1
        rhoC2 = 3.88968
    else if (x > 0.05 .and. r >= r0) then ! Omega2
        rhoC2 = 1
    else ! Omega3
        rhoC2 = 10
    end if
    
    end
    
    
    function u1C2(x,y)
    
    real x,y
    real u1C2
    
    if (x <= 0.05) then ! Omega1
        u1C2 = 11.2536
    else ! Omega2 & Omega3
        u1C2 = 0
    end if
    
    end
    
    
    function u2C2(x,y)
    
    real x,y
    real u2C2
    
    if (x <= 0.05) then ! Omega1
        u2C2 = 0
    else ! Omega2 & Omega3
        u2C2 = 0
    end if
    
    end
    
    
    function u3C2(x,y)
    
    real x,y
    real u3C2
    
    if (x <= 0.05) then ! Omega1
        u3C2 = 0
    else ! Omega2 & Omega3
        u3C2 = 0
    end if
    
    end
    
    
    function pC2(x,y)
    
    real x,y
    real pC2
    
    if (x <= 0.05) then ! Omega1
        pC2 = 167.345
    else ! Omega2 & Omega3
        pC2 = 1
    end if
    
    end
    
    
    function B1C2(x,y)
    
    real x,y
    real B1C2
    
    if (x <= 0.05) then ! Omega1
        B1C2 = 0
    else ! Omega2 & Omega3
        B1C2 = 0
    end if
    
    end
    
    
    
    function B2C2(x,y)
    
    real x,y
    real B2C2
    
    if (x <= 0.05) then ! Omega1
        B2C2 = 0.21826182
    else ! Omega2 & Omega3
        B2C2 = 0.56418958
    end if
    
    end
    
    
    function B3C2(x,y)
    
    real x,y
    real B3C2
    
    if (x <= 0.05) then ! Omega1
        B3C2 = -0.21826182
    else ! Omega2 & Omega3
        B3C2 = 0.56418958
    end if
    
    end