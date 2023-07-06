    function rhoC(x,y)
    
    real x,y,r,r0
    real rhoC
    
    r = ((x - 1.4)**2 + (y - 0.5)**2)**0.5
    r0 = 0.18
    
    if (x <= 1.2) then ! Omega1
        rhoC = 3.88968
    else if (x > 1.2 .and. r >= r0) then ! Omega2
        rhoC = 1
    else ! Omega3
        rhoC = 5
    end if
    
    end
    
    
    function u1C(x,y)
    
    real x,y
    real u1C
    
    if (x <= 1.2) then ! Omega1
        u1C = 0
    else ! Omega2 & Omega3
        u1C = -3.3156
    end if
    
    end
    
    
    function u2C(x,y)
    
    real x,y
    real u2C
    
    if (x <= 1.2) then ! Omega1
        u2C = 0
    else ! Omega2 & Omega3
        u2C = 0
    end if
    
    end
    
    
    function u3C(x,y)
    
    real x,y
    real u3C
    
    if (x <= 1.2) then ! Omega1
        u3C = -0.05234
    else ! Omega2 & Omega3
        u3C = 0
    end if
    
    end
    
    
    function pC(x,y)
    
    real x,y
    real pC
    
    if (x <= 1.2) then ! Omega1
        pC = 14.2641
    else ! Omega2 & Omega3
        pC = 0.04
    end if
    
    end
    
    
    function B1C(x,y)
    
    real x,y
    real B1C
    
    if (x <= 1.2) then ! Omega1
        B1C = 1
    else ! Omega2 & Omega3
        B1C = 1
    end if
    
    end
    
    
    
    function B2C(x,y)
    
    real x,y
    real B2C
    
    if (x <= 1.2) then ! Omega1
        B2C = 0
    else ! Omega2 & Omega3
        B2C = 0
    end if
    
    end
    
    
    function B3C(x,y)
    
    real x,y
    real B3C
    
    if (x <= 1.2) then ! Omega1
        B3C = 3.9353
    else ! Omega2 & Omega3
        B3C = 1
    end if
    
    end
    