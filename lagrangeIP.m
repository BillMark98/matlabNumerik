function [y] = lagrangeIP(t,ft,x)
    n = length(x);
    for index = 1 : n
        y(index) = lagrangeInterSum(t,ft,x(index));
    end
end

% evaluate the ith term of lagrange polynomial at data points t, with
% function value ft, at the point xi
function prod = lagrange_ith_term(i,t,ft,xi)
    n = length(t);
    prod = 1;
    for index = 1 : n
        if index ~= i
            prod = prod * (xi - t(index))/(t(i) - t(index));
        end
    end
    prod = prod * ft(i);
end

% evaluate the value of the interpolated polynomial at the given point xi
% given the data points specified in t, and data ft
function sum = lagrangeInterSum(t,ft,xi)
    n = length(t);
    sum = 0;
    for index = 1 : n
        sum = sum + lagrange_ith_term(index,t,ft,xi);
    end
end
