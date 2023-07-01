    subroutine HLL1(UR,UL,FR,FL,SR,SL,Uhat,Fhat)
    
    real UR(8),UL(8),FR(8),FL(8),SR,SL,Uhat(8),Fhat(8)
    
    if (SR < 0) then
        Fhat = FR
        Uhat = UR
    else if (SL > 0) then
        Fhat = FL
        Uhat = UL
    else
        Fhat = ( SR*FL - SL*FR + SL*SR*(UR - UL) )/(SR - SL)
        Uhat = ( SR*UR - SL*UL - (FR - FL) )/(SR - SL)
    end if
    
    end subroutine HLL1