    function pressure(rho,rhou,rhov,E,gamma)
    
    real rho,rhou,rhov,rhow,E,B1,B2,B3,gamma
    real pressure
    
    pressure = (gamma - 1)*(E - 0.5*(rhou**2 + rhov**2)/rho)
    
    end