function p = pressure(rho,rhou,E)

global gamma

p = (gamma - 1)*(E - 0.5*rhou^2/rho);

end