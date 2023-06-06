    subroutine HLL_Flux
    
    include 'com.txt'
    
    if (SR < 0) then
        Fhat1 = FR1
    else if (SL > 0) then
        Fhat1 = FL1
    else
        Fhat1 = ( SR*FL1 - SL*FR1 + SL*SR*(UR1 - UL1) )/(SR - SL)
    end if
    
    end subroutine HLL_Flux