function fu2 = f2(u1,u2,u3,u4,u5,u6)

global gamma
pstar = (gamma - 1)*(u4 - 0.5*(u2^2 + u3^2)/u1 - 0.5*(u5^2 + u6^2)) + 0.5*(u5^2 + u6^2);
fu2 = u2^2/u1 + pstar - u5^2;

end