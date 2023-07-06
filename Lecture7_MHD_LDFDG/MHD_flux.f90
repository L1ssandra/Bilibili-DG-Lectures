    subroutine MHD_flux(Uh,Fx,Fy)
    
    real Uh(8),Fx(8),Fy(8)
    real rho,u,v,w,p,E,B1,B2,B3,gamma,gamma1
    real S,T,K
      
    gamma = 5d0/3d0
    gamma1 = gamma - 1
    
    rho = Uh(1)
    u = Uh(2)/rho
    v = Uh(3)/rho
    w = Uh(4)/rho
    E = Uh(5)
    B1 = Uh(6)
    B2 = Uh(7)
    B3 = Uh(8)

    p = gamma1*(E - 0.5*rho*(u**2 + v**2 + w**2) - 0.5*(B1**2 + B2**2 + B3**2))
    S = p + 0.5*(B1**2 + B2**2 + B3**2)
    T = E + S
    K = u*B1 + v*B2 + w*B3
    
    Fx(1) = Uh(2)
    Fx(2) = rho*u**2 + S - B1**2
    Fx(3) = rho*u*v - B1*B2
    Fx(4) = rho*u*w - B1*B3
    Fx(5) = T*u - K*B1
    Fx(6) = 0
    Fx(7) = u*B2 - v*B1
    Fx(8) = u*B3 - w*B1
    
    Fy(1) = Uh(3)
    Fy(2) = rho*u*v - B1*B2
    Fy(3) = rho*v**2 + S - B2**2
    Fy(4) = rho*w*v - B3*B2
    Fy(5) = T*v - K*B2
    Fy(6) = v*B1 - u*B2
    Fy(7) = 0
    Fy(8) = v*B3 - w*B2
    
    end subroutine MHD_flux