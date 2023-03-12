% calculate_L2error_order.m

Table = zeros(4,2);
Nx0 = 20;

for n = 1:4
    Nx = Nx0*2^(n - 1);
    main
    Table(n,1) = L2_Error;
    if n > 1
        Table(n,2) = log(Table(n - 1,1)/Table(n,1))/log(2);
    end
end