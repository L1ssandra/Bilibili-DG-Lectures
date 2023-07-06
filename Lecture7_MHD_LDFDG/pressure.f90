    function pressure(rho,rhou,rhov,rhow,E,B1,B2,B3,gamma)
    
    real rho,rhou,rhov,rhow,E,B1,B2,B3,gamma
    real pressure
    
    pressure = (gamma - 1)*(E - 0.5*(rhou**2 + rhov**2 + rhow**2)/rho - 0.5*(B1**2 + B2**2 + B3**2))
    
    end