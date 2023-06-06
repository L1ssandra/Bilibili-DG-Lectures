    subroutine LF_Flux
    
    include 'com.txt'
    
    Fhat1 = 0.5d0*(FR1 + FL1 - max(abs(SR),abs(SL))*(UR1 - UL1))
    
    end subroutine LF_Flux