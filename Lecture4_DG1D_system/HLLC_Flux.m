function Fhat = HLLC_Flux(UR,UL,FR,FL,SR,SL)
global gamma
rhoL = UL(1); rhoR = UR(1);
rhouL = UL(2); rhouR = UR(2);
EL = UL(3); ER = UR(3);

if SL >= 0
    Fhat = FL;
elseif SR <= 0
    Fhat = FR;
else
    
    uL = rhouL/rhoL;
    uR = rhouR/rhoR;
    
    pL = (gamma - 1)*(EL - 0.5*rhoL*uL^2);
    pR = (gamma - 1)*(ER - 0.5*rhoR*uR^2);
    
    Sstar = (pR - pL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR))/(rhoL*(SL - uL) - rhoR*(SR - uR));
    
    Dstar = [0,1,Sstar];
    
    if SL < 0 && Sstar >= 0
        Fhat = ( Sstar*(SL*UL - FL) + SL*(pL + rhoL*(SL - uL)*(Sstar - uL))*Dstar )/(SL - Sstar);
    else
        Fhat = ( Sstar*(SR*UR - FR) + SR*(pR + rhoR*(SR - uR)*(Sstar - uR))*Dstar )/(SR - Sstar);
    end
end

end