    subroutine HLLC_Flux
    
    include 'com.txt'
    
    real Sstar,rhoR,rhoL,uRbot,uLbot,PtotR,PtotL,ER,EL,u1R,u1L,u2R,u2L,u3R,u3L,B1R,B1L,B2R,B2L,B3R,B3L
    real rhoRstar,rhoLstar,B1star,B2star,B3star,u1star,u2star,u1Rstar,u1Lstar,u2Rstar,u2Lstar,u3Rstar,u3Lstar
    real BRbot,BLbot,Bstarbot,Ptotstar,ERstar,ELstar
    
    if (SR < 0) then
        Fhat1 = FR1
        Ustar = UR1
    else if (SL > 0) then
        Fhat1 = FL1
        Ustar = UL1
    else
        Ustar = ( SR*UR1 - SL*UL1 + FL1 - FR1 )/(SR - SL)
        
        rhoR = UR1(1)
        u1R = UR1(2)/rhoR
        u2R = UR1(3)/rhoR
        u3R = UR1(4)/rhoR
        ER = UR1(5)
        B1R = UR1(6)
        B2R = UR1(7)
        B3R = UR1(8)
        
        rhoL = UL1(1)
        u1L = UL1(2)/rhoL
        u2L = UL1(3)/rhoL
        u3L = UL1(4)/rhoL
        EL = UL1(5)
        B1L = UL1(6)
        B2L = UL1(7)
        B3L = UL1(8)
        
        B1star = Ustar(6)
        B2star = Ustar(7)
        B3star = Ustar(8)
        
        if (direction == 1) then
            Sstar = Ustar(2)/Ustar(1)
            uRbot = u1R
            uLbot = u1L
            BRbot = B1R
            BLbot = B1L
            Bstarbot = B1star
        else if (direction == 2) then
            Sstar = Ustar(3)/Ustar(1)
            uRbot = u2R
            uLbot = u2L
            BRbot = B2R
            BLbot = B2L
            Bstarbot = B2star
        end if
        
        rhoRstar = rhoR*(SR - uRbot)/(SR - Sstar)
        rhoLstar = rhoL*(SL - uLbot)/(SL - Sstar)
        
        PtotR = gamma1*(ER - 0.5d0*rhoR*(u1R**2 + u2R**2 + u3R**2))
        PtotL = gamma1*(EL - 0.5d0*rhoL*(u1L**2 + u2L**2 + u3L**2))
        
        PtotRstar = PtotR + rhoR*(SR - uRbot)*(Sstar - uRbot) + Bstarbot**2 - BRbot**2
        PtotLstar = PtotL + rhoL*(SL - uLbot)*(Sstar - uLbot) + Bstarbot**2 - BLbot**2
        
        if (direction == 1) then
            URstar(2) = (UR1(2)*(SR - uRbot) + ptotRstar - ptotR + BRbot*B1R - Bstarbot*B1star)/(SR - Sstar)
            URstar(3) = (UR1(3)*(SR - uRbot)                     + BRbot*B2R - Bstarbot*B2star)/(SR - Sstar)
            ULstar(2) = (UL1(2)*(SL - uLbot) + ptotLstar - ptotL + BLbot*B1L - Bstarbot*B1star)/(SL - Sstar)
            ULstar(3) = (UL1(3)*(SL - uLbot)                     + BLbot*B2L - Bstarbot*B2star)/(SL - Sstar)
        else if (direction == 2) then
            URstar(2) = (UR1(2)*(SR - uRbot)                     + BRbot*B1R - Bstarbot*B1star)/(SR - Sstar)
            URstar(3) = (UR1(3)*(SR - uRbot) + ptotRstar - ptotR + BRbot*B2R - Bstarbot*B2star)/(SR - Sstar)
            ULstar(2) = (UL1(2)*(SL - uLbot)                     + BLbot*B1L - Bstarbot*B1star)/(SL - Sstar)
            ULstar(3) = (UL1(3)*(SL - uLbot) + ptotLstar - ptotL + BLbot*B2L - Bstarbot*B2star)/(SL - Sstar)
        end if
        URstar(4) = (UR1(4)*(SR - uRbot) + BRbot*B3R - Bstarbot*B3star)/(SR - Sstar)
        ULstar(4) = (UL1(4)*(SL - uLbot) + BLbot*B3L - Bstarbot*B3star)/(SL - Sstar)
        
        u1Rstar = URstar(2)/rhoRstar
        u2Rstar = URstar(3)/rhoRstar
        u3Rstar = URstar(4)/rhoRstar
        
        u1Lstar = ULstar(2)/rhoLstar
        u2Lstar = ULstar(3)/rhoLstar
        u3Lstar = ULstar(4)/rhoLstar
        
        ERstar = ( (SR - uRbot)*ER - PtotR*uRbot + PtotRstar*Sstar + BRbot*(u1R*B1R + u2R*B2R + u3R*B3R) - Bstarbot*(u1Rstar*B1star + u2Rstar*B2star + u3Rstar*B3star) )/(SR - Sstar)
        ELstar = ( (SL - uLbot)*EL - PtotL*uLbot + PtotLstar*Sstar + BLbot*(u1L*B1L + u2L*B2L + u3L*B3L) - Bstarbot*(u1Lstar*B1star + u2Lstar*B2star + u3Lstar*B3star) )/(SL - Sstar)
        
        URstar(1) = rhoRstar
        URstar(5) = ERstar
        URstar(6:8) = Ustar(6:8)
        
        ULstar(1) = rhoLstar
        ULstar(5) = ELstar
        ULstar(6:8) = Ustar(6:8)
    
        if (Sstar > 0) then
            Fhat1 = FL1 + SL*(ULstar - UL1)
            Ustar = ULstar
        else
            Fhat1 = FR1 + SR*(URstar - UR1)
            Ustar = URstar
        end if
        
    end if
    
    end subroutine HLLC_Flux