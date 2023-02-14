% L2Pro3.m

N0 = 10;
Table = zeros(4,2);
for m = 1:4
    N = N0*2^(m - 1);
    L2Pro2
    Table(m,1) = L2_Error;
    if m > 1
        Table(m,2) = log(Table(m - 1,1)/Table(m,1))/log(2);
    end
end