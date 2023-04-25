function [SR,SL] = wavespeed(u)

global gamma

p = (gamma - 1)*(u(3) - 0.5*u(2)^2/u(1));
c = sqrt(gamma*p/u(1));

SR = u(2)/u(1) + c;
SL = u(2)/u(1) - c;

end