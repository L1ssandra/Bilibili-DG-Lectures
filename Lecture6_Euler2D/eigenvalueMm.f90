    subroutine eigenvalueMm(Amax,Amin,rho,rhou,rhov,E,n1,n2)
    
    include 'com.txt'
    real u,v,w,p,c,BP,Bn,un
    
    u = rhou/rho
    v = rhov/rho
    
    un = u*n1 + v*n2
    
    p = gamma1*(E - 0.5d0*rho*(u**2 + v**2))
    
    c = sqrt(abs(gamma*p/rho))
    
    Amax = un + c
    Amin = un - c
    
    end subroutine eigenvalueMm
    
    