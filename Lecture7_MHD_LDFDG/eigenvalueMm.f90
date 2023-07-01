    subroutine eigenvalueMm(Amax,Amin,rho,rhou,rhov,rhow,E,B1,B2,B3,n1,n2)
    
    include 'com.txt'
    real u,v,w,p,c,BP,Bn,un,cf,n3
    
    n3 = 0
    
    u = rhou/rho
    v = rhov/rho
    w = rhow/rho
    
    BP = B1**2 + B2**2 + B3**2
    Bn = B1*n1 + B2*n2 + B3*n3
    un = u*n1 + v*n2 + w*n3
    
    p = gamma1*(E - 0.5d0*rho*(u**2 + v**2 + w**2) - 0.5d0*BP)
    
    c = sqrt(abs(gamma*p/rho))
    
    cf = sqrt(abs( 0.5d0*(c**2 + BP/rho + sqrt((c**2 + BP/rho)**2 - 4*c**2*Bn**2/rho) ) ))
    
    Amax = un + cf
    Amin = un - cf
    
    end subroutine eigenvalueMm
    
    