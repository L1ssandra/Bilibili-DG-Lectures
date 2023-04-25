function fhat = HLL_Flux(uR,uL,fR,fL,SR,SL)

if SL >= 0
    fhat = fL;
elseif SR <= 0
    fhat = fR;
else
    fhat = (SR*fL - SL*fR + SL*SR*(uR - uL))/(SR - SL);
end

end