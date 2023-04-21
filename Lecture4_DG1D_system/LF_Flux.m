function fhat = LF_Flux(uR,uL,fR,fL,SR,SL)

alpha = max(abs(SR),abs(SL));
fhat = 0.5*(fR + fL) - 0.5*alpha*(uR - uL);

end