% calculate_L2error_order.m
clear;clc
Table = zeros(4,6);
Nx0 = 5;

global Nx

for ii = 1:4
    Nx = Nx0*2^(ii - 1);
    main
    Table(ii,1) = L2_Error(1);
    Table(ii,3) = L2_Error(2);
    Table(ii,5) = L2_Error(3);
    if ii > 1
        Table(ii,2) = log(Table(ii - 1,1)/Table(ii,1))/log(2);
        Table(ii,4) = log(Table(ii - 1,3)/Table(ii,3))/log(2);
        Table(ii,6) = log(Table(ii - 1,5)/Table(ii,5))/log(2);
    end
end