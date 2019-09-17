% calculate the sum of the absolute value all the Lagrange poly of degree n
% for data points ti at a specific point t
function sum = lagrangeInter(n,ti,t)
    sum = 0;
    for index = 0 : n
        sum  = sum + abs(lagrange_ith_Inter(index,n,ti,t));
    end
end