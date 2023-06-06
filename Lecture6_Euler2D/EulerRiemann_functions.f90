    function rho1(x,y)
    
    real x,y
    real rho1
    
    if (x > 0.8 .and. y > 0.8) then
        rho1 = 1.5
    else if (x <= 0.8 .and. y > 0.8) then
        rho1 = 0.5323
    else if (x <= 0.8 .and. y <= 0.8) then
        rho1 = 0.138
    else
        rho1 = 0.5323
    end if
    
    end
    !*********************************************
    function ux1(x,y)
    
    real x,y
    real ux1
    
    if (x > 0.8 .and. y > 0.8) then
        ux1 = 0
    else if (x <= 0.8 .and. y > 0.8) then
        ux1 = 1.206
    else if (x <= 0.8 .and. y <= 0.8) then
        ux1 = 1.206
    else
        ux1 = 0
    end if
    
    end
    !*********************************************
    function uy1(x,y)
    
    real x,y
    real uy1
    
    if (x > 0.8 .and. y > 0.8) then
        uy1 = 0
    else if (x <= 0.8 .and. y > 0.8) then
        uy1 = 0
    else if (x <= 0.8 .and. y <= 0.8) then
        uy1 = 1.206
    else
        uy1 = 1.206
    end if
    
    end
    !*********************************************
    function pp1(x,y)
    
    real x,y
    real pp1
    
    if (x > 0.8 .and. y > 0.8) then
        pp1 = 1.5
    else if (x <= 0.8 .and. y > 0.8) then
        pp1 = 0.3
    else if (x <= 0.8 .and. y <= 0.8) then
        pp1 = 0.029
    else
        pp1 = 0.3
    end if
    
    end
