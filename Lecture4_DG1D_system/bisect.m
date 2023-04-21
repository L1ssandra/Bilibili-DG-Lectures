function xstar = bisect(f,a,b,tol)

if f(a)*f(b) > tol
    fprintf('%d, %d\n',f(a),f(b))
    %error('f(a)f(b) < 0 not satisfied!')
elseif abs(f(a)) < tol
    b = a;
elseif abs(f(b)) < tol
    a = b;
end

fa = f(a);

while (b - a)/2 >= tol
    c = (a + b)/2;
    fc = f(c);
    if fc == 0
        break
    end
    if sign(fc) * sign(fa) < 0
        b = c;
    else
        a = c;
        fa = fc;
    end
end

xstar = a;

end