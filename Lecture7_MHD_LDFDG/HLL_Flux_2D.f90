    subroutine HLL_Flux_2D
    
    include 'com.txt'
    
    real alphaxRU,alphaxLU,alphaxRD,alphaxLD
    real alphayRU,alphayLU,alphayRD,alphayLD
    real alphax2D,alphay2D
    real EzRU,EzLU,EzRD,EzLD
    real BxRU,BxLU,BxRD,BxLD
    real ByRU,ByLU,ByRD,ByLD
    real EzR1,EzL1,EzU1,EzD1
    real EzRstar,EzLstar,EzUstar,EzDstar,EzStarStar
    real Ustarstar(NumEq)
    real FRstar(NumEq),FLstar(NumEq),FUstar(NumEq),FDstar(NumEq)
    real FxRU(NumEq),FyRU(NumEq),FxLU(NumEq),FyLU(NumEq),FxRD(NumEq),FyRD(NumEq),FxLD(NumEq),FyLD(NumEq)
    
    call eigenvalueMm(alphaxRU,betaxRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),1,0)
    call eigenvalueMm(alphayRU,betayRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),0,1)
    call eigenvalueMm(alphaxLU,betaxLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),1,0)
    call eigenvalueMm(alphayLU,betayLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),0,1)
    call eigenvalueMm(alphaxRD,betaxRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),1,0)
    call eigenvalueMm(alphayRD,betayRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),0,1)
    call eigenvalueMm(alphaxLD,betaxLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),1,0)
    call eigenvalueMm(alphayLD,betayLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),0,1)
    
    SR = max(alphaxRU,alphaxLU,alphaxRD,alphaxLD)
    SL = min(betaxRU,betaxLU,betaxRD,betaxLD)
    SU = max(alphayRU,alphayLU,alphayRD,alphayLD)
    SD = min(betayRU,betayLU,betayRD,betayLD)
    
    call MHD_flux(URU1,FxRU,FyRU)
    call MHD_flux(ULU1,FxLU,FyLU)
    call MHD_flux(URD1,FxRD,FyRD)
    call MHD_flux(ULD1,FxLD,FyLD)
    
    BxRU = URU1(6)
    BxLU = ULU1(6)
    BxRD = URD1(6)
    BxLD = ULD1(6)
    
    ByRU = URU1(7)
    ByLU = ULU1(7)
    ByRD = URD1(7)
    ByLD = ULD1(7)
    
    call HLL1(URU1,ULU1,FxRU,FxLU,SR,SL,UUstar,FUstar)
    call HLL1(URD1,ULD1,FxRD,FxLD,SR,SL,UDstar,FDstar)
    call HLL1(URU1,URD1,FyRU,FyRD,SU,SD,URstar,FRstar)
    call HLL1(ULU1,ULD1,FyLU,FyLD,SU,SD,ULstar,FLstar)
    
    BxUstar = UUstar(6)
    BxDstar = UDstar(6)
    ByRstar = URstar(7)
    ByLstar = ULstar(7)
    
    EzRstar = FRstar(6)
    EzLstar = FLstar(6)
    EzUstar = -FUstar(7)
    EzDstar = -FDstar(7)
    
    BxStarStar = ( 2*SR*SU*BxRU - 2*SL*SU*BxLU + 2*SL*SD*BxLD - 2*SR*SD*BxRD - SR*(EzRU - EzRD) + SL*(EzLU - EzLD) - (SR - SL)*(EzUstar - EzDstar) )/( 2*(SR - SL)*(SU - SD) )
    ByStarStar = ( 2*SR*SU*ByRU - 2*SL*SU*ByLU + 2*SL*SD*ByLD - 2*SR*SD*ByRD + SU*(EzRU - EzLU) - SD*(EzRD - EzLD) + (SU - SD)*(EzRstar - EzLstar) )/( 2*(SR - SL)*(SU - SD) )
    
    EzStarStar = 0.25*(EzRstar + EzLstar + EzUstar + EzDstar) - 0.25*SU*(BxUstar - BxStarStar) - 0.25*SD*(BxDstar - BxStarStar) + 0.25*SR*(ByRstar - ByStarStar) + 0.25*SL*(ByLstar - ByStarStar)
    
    !print *," Ez "
    !print *,URU1
    !print *," "
    !print *,ULU1
    !print *," "
    !print *,URD1
    !print *," "
    !print *,ULD1
    !print *," "
    
    !EzR1 = 0.5*(EzRU + EzRD - alphay2D*(B1RU - B1RD) )
    !EzL1 = 0.5*(EzLU + EzLD - alphay2D*(B1LU - B1LD) )
    !EzU1 = 0.5*(EzRU + EzLU + alphax2D*(B2RU - B2LU) )
    !EzD1 = 0.5*(EzRD + EzLD + alphax2D*(B2RD - B2LD) )
    
    !print *,EzRstar
    !print *,EzLstar
    !print *,EzUstar
    !print *,EzDstar
    !print *," "
    !print *,EzUD(21,20,1)
    !print *,EzUD(20,20,NumGLP)
    !print *,EzRL(20,21,1)
    !print *,EzRL(20,20,NumGLP)
    !print *," "
    
    !print *,EzRU,EzLU,EzRD,EzLD
    
    !print *,(0.5*(B1RU + B1LU) - 0.5*(B1RD + B1LD)),alphax2D*(0.5*(B2RU + B2RD) - 0.5*(B2LU + B2LD))
    !print *," "
    !print *,Ezhat,alphay2D,alphax2D
    !print *,0.25*(EzRL(20,20,NumGLP) + EzRL(20,21,1) + EzUD(20,20,NumGLP) + EzUD(21,20,1))
    
    if (SL > 0) then
        Ezhat = EzLstar
    else if (SR < 0) then
        Ezhat = EzRstar
    else if (SD > 0) then
        Ezhat = EzDstar
    else if (SU < 0) then
        Ezhat = EzUstar
    else
        Ezhat = EzStarStar
    end if
    
    end subroutine HLL_Flux_2D