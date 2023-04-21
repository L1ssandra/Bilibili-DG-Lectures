function fu2 = f2(u1,u2,u3)

global gamma
p = (gamma - 1)*(u3 - 0.5*u2^2/u1);
fu2 = u2^2/u1 + p;

end