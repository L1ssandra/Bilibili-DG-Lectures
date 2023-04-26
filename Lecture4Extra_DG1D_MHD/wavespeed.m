function [SR,SL] = wavespeed(u)

global gamma

p = (gamma - 1)*(u(4) - 0.5*(u(2)^2 + u(3)^2)/u(1) - 0.5*(u(5)^2 + u(6)^2));
c = sqrt(gamma*p/u(1));
BP = u(5)^2 + u(6)^2;

cf = sqrt(abs( 0.5*(c^2 + BP/u(1) + sqrt((c^2 + BP/u(1))^2 - 4*c^2*u(5)^2/u(1)) ) ));

SR = u(2)/u(1) + cf;
SL = u(2)/u(1) - cf;

end