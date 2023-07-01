    subroutine HLLD_Flux
    
    include 'com.txt'
    
    real Sstar,rhoR,rhoL,uRbot,uLbot,PtotR,PtotL,ER,EL,u1R,u1L,u2R,u2L,u3R,u3L,B1R,B1L,B2R,B2L,B3R,B3L
    real rhoRstar,rhoLstar,B1star,B2star,B3star,u1star,u2star,u1Rstar,u1Lstar,u2Rstar,u2Lstar,u3Rstar,u3Lstar
    real BRbot,BLbot,Bstarbot,ptotRstar,ptotLstar,ERstar,ELstar,SRstar,SLstar
    real B1Lstar,B1Rstar,B2Lstar,B2Rstar,B3Lstar,B3Rstar
    real AL,AR,CL,CR,GL,GR,DL,DR,rhoRstarstar,rhoLstarstar
    real u1starstar,u2starstar,u3starstar
    real B1starstar,B2starstar,B3starstar
    real ERstarstar,ELstarstar
    
    if (SR < 0) then
        Fhat1 = FR1
        Ustar = UR1
    else if (SL > 0) then
        Fhat1 = FL1
        Ustar = UL1
    else
        ! Stage 1
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
        
        if (direction == 1) then
            Sstar = Ustar(2)/Ustar(1)
            uRbot = u1R
            uLbot = u1L
            BRbot = B1R
            BLbot = B1L
            Bstarbot = Ustar(6)
        else if (direction == 2) then
            Sstar = Ustar(3)/Ustar(1)
            uRbot = u2R
            uLbot = u2L
            BRbot = B2R
            BLbot = B2L
            Bstarbot = Ustar(7)
        end if
        
        rhoRstar = rhoR*(SR - uRbot)/(SR - Sstar)
        rhoLstar = rhoL*(SL - uLbot)/(SL - Sstar)
        
        SRstar = Sstar + (Bstarbot**2/rhoRstar)**0.5
        SLstar = Sstar - (Bstarbot**2/rhoLstar)**0.5
        
        PtotR = gamma1*(ER - 0.5d0*rhoR*(u1R**2 + u2R**2 + u3R**2))
        PtotL = gamma1*(EL - 0.5d0*rhoL*(u1L**2 + u2L**2 + u3L**2))
        
        PtotRstar = PtotR + rhoR*(SR - uRbot)*(Sstar - uRbot) + Bstarbot**2 - BRbot**2
        PtotLstar = PtotL + rhoL*(SL - uLbot)*(Sstar - uLbot) + Bstarbot**2 - BLbot**2
        
        AR = rhoR*(SR - uRbot)**2 - Bstarbot*BRbot
        AL = rhoL*(SL - uLbot)**2 - Bstarbot*BLbot
        
        CR = rhoR*(SR - uRbot)*(BRbot - Bstarbot)
        CL = rhoL*(SL - uLbot)*(BLbot - Bstarbot)
        
        DR = rhoR*(SR - uRbot)*(SR - Sstar) - Bstarbot**2 + 1e-15
        DL = rhoL*(SL - uLbot)*(SL - Sstar) - Bstarbot**2 + 1e-15
        
        if (direction == 1) then
            B1Rstar = Bstarbot
            B1Lstar = Bstarbot
            B2Rstar = (AR*B2R + CR*u2R)/DR
            B2Lstar = (AL*B2L + CL*u2L)/DL
        else if (direction == 2) then
            B1Rstar = (AR*B1R + CR*u1R)/DR
            B1Lstar = (AL*B1L + CL*u1L)/DL
            B2Rstar = Bstarbot
            B2Lstar = Bstarbot
        end if
        B3Rstar = (AR*B3R + CR*u3R)/DR
        B3Lstar = (AL*B3L + CL*u3L)/DL
        
        if (direction == 1) then
            URstar(2) = (UR1(2)*(SR - uRbot) + ptotRstar - ptotR + BRbot*B1R - Bstarbot*B1Rstar)/(SR - Sstar)
            URstar(3) = (UR1(3)*(SR - uRbot)                     + BRbot*B2R - Bstarbot*B2Rstar)/(SR - Sstar)
            ULstar(2) = (UL1(2)*(SL - uLbot) + ptotLstar - ptotL + BLbot*B1L - Bstarbot*B1Lstar)/(SL - Sstar)
            ULstar(3) = (UL1(3)*(SL - uLbot)                     + BLbot*B2L - Bstarbot*B2Lstar)/(SL - Sstar)
        else if (direction == 2) then
            URstar(2) = (UR1(2)*(SR - uRbot)                     + BRbot*B1R - Bstarbot*B1Rstar)/(SR - Sstar)
            URstar(3) = (UR1(3)*(SR - uRbot) + ptotRstar - ptotR + BRbot*B2R - Bstarbot*B2Rstar)/(SR - Sstar)
            ULstar(2) = (UL1(2)*(SL - uLbot)                     + BLbot*B1L - Bstarbot*B1Lstar)/(SL - Sstar)
            ULstar(3) = (UL1(3)*(SL - uLbot) + ptotLstar - ptotL + BLbot*B2L - Bstarbot*B2Lstar)/(SL - Sstar)
        end if
        URstar(4) = (UR1(4)*(SR - uRbot) + BRbot*B3R - Bstarbot*B3Rstar)/(SR - Sstar)
        ULstar(4) = (UL1(4)*(SL - uLbot) + BLbot*B3L - Bstarbot*B3Lstar)/(SL - Sstar)
        
        u1Rstar = URstar(2)/rhoRstar
        u2Rstar = URstar(3)/rhoRstar
        u3Rstar = URstar(4)/rhoRstar
        
        u1Lstar = ULstar(2)/rhoLstar
        u2Lstar = ULstar(3)/rhoLstar
        u3Lstar = ULstar(4)/rhoLstar
        
        ERstar = ( (SR - uRbot)*ER - PtotR*uRbot + PtotRstar*Sstar + BRbot*(u1R*B1R + u2R*B2R + u3R*B3R) - Bstarbot*(u1Rstar*B1Rstar + u2Rstar*B2Rstar + u3Rstar*B3Rstar) )/(SR - Sstar)
        ELstar = ( (SL - uLbot)*EL - PtotL*uLbot + PtotLstar*Sstar + BLbot*(u1L*B1L + u2L*B2L + u3L*B3L) - Bstarbot*(u1Lstar*B1Lstar + u2Lstar*B2Lstar + u3Lstar*B3Lstar) )/(SL - Sstar)
        
        URstar(1) = rhoRstar
        URstar(5) = ERstar
        URstar(6:8) = Ustar(6:8)
        
        ULstar(1) = rhoLstar
        ULstar(5) = ELstar
        ULstar(6:8) = Ustar(6:8)
        
        if ((abs(Bstarbot) > 1e-4) .and. (SRstar - SLstar > 1e-4)) then
            if (SLstar > 0) then
                Fhat1 = FL1 + SL*(ULstar - UL1)
                Ustar = ULstar
            else if (SRstar < 0) then
                Fhat1 = FR1 + SR*(URstar - UR1)
                Ustar = URstar
            else
                ! Stage 2
                FR1 = FR1 + SR*(URstar - UR1)
                FL1 = FL1 + SL*(ULstar - UL1)
            
                rhoRstarstar = rhoRstar
                rhoLstarstar = rhoLstar
            
                u1starstar = (abs(rhoLstar)**0.5*u1Lstar + abs(rhoRstar)**0.5*u1Rstar + sign(1d0,Bstarbot)*(B1Rstar - B1Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
                u2starstar = (abs(rhoLstar)**0.5*u2Lstar + abs(rhoRstar)**0.5*u2Rstar + sign(1d0,Bstarbot)*(B2Rstar - B2Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
                u3starstar = (abs(rhoLstar)**0.5*u3Lstar + abs(rhoRstar)**0.5*u3Rstar + sign(1d0,Bstarbot)*(B3Rstar - B3Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
            
                B1starstar = (abs(rhoLstar)**0.5*B1Lstar + abs(rhoRstar)**0.5*B1Rstar + sign(1d0,Bstarbot)*abs(rhoLstar*rhoRstar)**0.5*(u1Rstar - u1Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
                B2starstar = (abs(rhoLstar)**0.5*B2Lstar + abs(rhoRstar)**0.5*B2Rstar + sign(1d0,Bstarbot)*abs(rhoLstar*rhoRstar)**0.5*(u2Rstar - u2Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
                B3starstar = (abs(rhoLstar)**0.5*B3Lstar + abs(rhoRstar)**0.5*B3Rstar + sign(1d0,Bstarbot)*abs(rhoLstar*rhoRstar)**0.5*(u3Rstar - u3Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
            
                ERstarstar = ERstar + abs(rhoRstar)**0.5*(u1Rstar*B1Rstar + u2Rstar*B2Rstar + u3Rstar*B3Rstar - u1starstar*B1starstar - u2starstar*B2starstar - u3starstar*B3starstar)*sign(1d0,Bstarbot)
                ELstarstar = ELstar - abs(rhoLstar)**0.5*(u1Lstar*B1Lstar + u2Lstar*B2Lstar + u3Lstar*B3Lstar - u1starstar*B1starstar - u2starstar*B2starstar - u3starstar*B3starstar)*sign(1d0,Bstarbot)
                
                URstarstar(1) = rhoRstarstar
                URstarstar(2) = rhoRstarstar*u1starstar
                URstarstar(3) = rhoRstarstar*u2starstar
                URstarstar(4) = rhoRstarstar*u3starstar
                URstarstar(5) = ERstarstar
                URstarstar(6) = B1starstar
                URstarstar(7) = B2starstar
                URstarstar(8) = B3starstar
            
                ULstarstar(1) = rhoLstarstar
                ULstarstar(2) = rhoLstarstar*u1starstar
                ULstarstar(3) = rhoLstarstar*u2starstar
                ULstarstar(4) = rhoLstarstar*u3starstar
                ULstarstar(5) = ELstarstar
                ULstarstar(6) = B1starstar
                ULstarstar(7) = B2starstar
                ULstarstar(8) = B3starstar
            
                if (Sstar > 0) then
                    Fhat1 = FL1 + SLstar*(ULstarstar - ULstar)
                    Ustar = ULstarstar
                else
                    Fhat1 = FR1 + SRstar*(URstarstar - URstar)
                    Ustar = URstarstar
                end if
                
            end if
        else
            if (Sstar > 0) then
                Fhat1 = FL1 + SL*(ULstar - UL1)
                Ustar = ULstar
            else
                Fhat1 = FR1 + SR*(URstar - UR1)
                Ustar = URstar
            end if
        end if
        
    end if
    
    end subroutine HLLD_Flux