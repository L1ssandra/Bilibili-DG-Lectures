function fu3 = f3(u1,u2,u3)

global gamma
p = (gamma - 1)*(u3 - 0.5*u2^2/u1);
fu3 = u2/u1*(u3 + p);

end