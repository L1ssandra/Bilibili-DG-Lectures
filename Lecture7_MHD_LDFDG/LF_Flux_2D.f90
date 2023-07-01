    subroutine LF_Flux_2D
    
    include 'com.txt'
    
    real alphaxRU,alphaxLU,alphaxRD,alphaxLD
    real alphayRU,alphayLU,alphayRD,alphayLD
    real alphax2D,alphay2D
    real EzRU,EzLU,EzRD,EzLD
    real B1RU,B1LU,B1RD,B1LD
    real B2RU,B2LU,B2RD,B2LD
    real EzR1,EzL1,EzU1,EzD1
    
    !alphax2D = max(alphaxRU,alphaxLU,alphaxRD,alphaxLD)
    !alphay2D = max(alphayRU,alphayLU,alphayRD,alphayLD)
    alphax2D = max(abs(URU1(2)/URU1(1)),abs(ULU1(2)/ULU1(1)),abs(URD1(2)/URD1(1)),abs(ULD1(2)/ULD1(1)))
    alphay2D = max(abs(URU1(3)/URU1(1)),abs(ULU1(3)/ULU1(1)),abs(URD1(3)/URD1(1)),abs(ULD1(3)/ULD1(1)))
    
    call calculate_Ez(EzRU,URU1(2)/URU1(1),URU1(3)/URU1(1),URU1(6),URU1(7))
    call calculate_Ez(EzLU,ULU1(2)/ULU1(1),ULU1(3)/ULU1(1),ULU1(6),ULU1(7))
    call calculate_Ez(EzRD,URD1(2)/URD1(1),URD1(3)/URD1(1),URD1(6),URD1(7))
    call calculate_Ez(EzLD,ULD1(2)/ULD1(1),ULD1(3)/ULD1(1),ULD1(6),ULD1(7))
    
    !print *," Ez "
    !print *,URU1
    !print *," "
    !print *,ULU1
    !print *," "
    !print *,URD1
    !print *," "
    !print *,ULD1
    !print *," "
    
    B1RU = URU1(6)
    B1LU = ULU1(6)
    B1RD = URD1(6)
    B1LD = ULD1(6)
    
    B2RU = URU1(7)
    B2LU = ULU1(7)
    B2RD = URD1(7)
    B2LD = ULD1(7)
    
    !EzR1 = 0.5*(EzRU + EzRD - alphay2D*(B1RU - B1RD) )
    !EzL1 = 0.5*(EzLU + EzLD - alphay2D*(B1LU - B1LD) )
    !EzU1 = 0.5*(EzRU + EzLU + alphax2D*(B2RU - B2LU) )
    !EzD1 = 0.5*(EzRD + EzLD + alphax2D*(B2RD - B2LD) )
    
    !print *,EzR1
    !print *,EzL1
    !print *,EzU1
    !print *,EzD1
    !print *," "
    !print *,EzUD(21,20,1)
    !print *,EzUD(20,20,NumGLP)
    !print *,EzRL(20,21,1)
    !print *,EzRL(20,20,NumGLP)
    !print *," "
    
    !print *,EzRU,EzLU,EzRD,EzLD
    Ezhat = 0.25*(EzRU + EzLU + EzRD + EzLD) - 0.25*alphay2D*(0.5*(B1RU + B1LU) - 0.5*(B1RD + B1LD)) + 0.25*alphax2D*(0.5*(B2RU + B2RD) - 0.5*(B2LU + B2LD))
    
    !print *,(0.5*(B1RU + B1LU) - 0.5*(B1RD + B1LD)),alphax2D*(0.5*(B2RU + B2RD) - 0.5*(B2LU + B2LD))
    !print *," "
    !print *,Ezhat,alphay2D,alphax2D
    !print *,0.25*(EzRL(20,20,NumGLP) + EzRL(20,21,1) + EzUD(20,20,NumGLP) + EzUD(21,20,1))
    
    end subroutine LF_Flux_2D