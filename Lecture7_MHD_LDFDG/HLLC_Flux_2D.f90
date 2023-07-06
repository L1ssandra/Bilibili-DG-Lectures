    subroutine HLLC_Flux_2D
    
    include 'com.txt'
    
    real alphaxRU,alphaxLU,alphaxRD,alphaxLD
    real alphayRU,alphayLU,alphayRD,alphayLD
    real alphax2D,alphay2D,SR1,SL1
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
    
    direction = 1
    
    UR1 = URU1
    UL1 = ULU1
    FR1 = FxRU
    FL1 = FxLU
    call HLLC_Flux
    FUstar = Fhat1
    UUstar = Ustar
    
    UR1 = URD1
    UL1 = ULD1
    FR1 = FxRD
    FL1 = FxLD
    call HLLC_Flux
    FDstar = Fhat1
    UDstar = Ustar
    
    SR1 = SR
    SL1 = SL
    SR = SU
    SL = SD
    direction = 2
    
    UR1 = URU1
    UL1 = URD1
    FR1 = FyRU
    FL1 = FyRD
    call HLLC_Flux
    FRstar = Fhat1
    URstar = Ustar
    
    UR1 = ULU1
    UL1 = ULD1
    FR1 = FyLU
    FL1 = FyLD
    call HLLC_Flux
    FLstar = Fhat1
    ULstar = Ustar
    
    SR = SR1
    SL = SL1
    
    BxUstar = UUstar(6)
    BxDstar = UDstar(6)
    ByRstar = URstar(7)
    ByLstar = ULstar(7)
    
    EzRstar = FRstar(6)
    EzLstar = FLstar(6)
    EzUstar = -FUstar(7)
    EzDstar = -FDstar(7)
    
    !BxStarStar = ( 2*SR*SU*BxRU - 2*SL*SU*BxLU + 2*SL*SD*BxLD - 2*SR*SD*BxRD - SR*(EzRU - EzRD) + SL*(EzLU - EzLD) - (SR - SL)*(EzUstar - EzDstar) )/( 2*(SR - SL)*(SU - SD) )
    !ByStarStar = ( 2*SR*SU*ByRU - 2*SL*SU*ByLU + 2*SL*SD*ByLD - 2*SR*SD*ByRD + SU*(EzRU - EzLU) - SD*(EzRD - EzLD) + (SU - SD)*(EzRstar - EzLstar) )/( 2*(SR - SL)*(SU - SD) )
    
    EzStarStar = 0.25*(EzRstar + EzLstar + EzUstar + EzDstar) !- 0.25*SU*(BxUstar - BxStarStar) - 0.25*SD*(BxDstar - BxStarStar) + 0.25*SR*(ByRstar - ByStarStar) + 0.25*SL*(ByLstar - ByStarStar)
    
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
    
    end subroutine HLLC_Flux_2D