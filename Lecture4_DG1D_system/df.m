function dfu = df(u)

dfu = (-40.*u + 8).*u.^2./(4.*u.^2 + (1 - u).^2).^2 + 8.*u./(4*u.^2 + (1 - u).^2);

end