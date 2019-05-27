function [y] = lagrangePolynomial(i,t,x)
% given a vector of points x
% calculate the i-th lagrange polynomial (begins at 0)
% with data points t
y = ones(size(x));
n = length(t);

for j = 1 : n
    if((i+1) ~= j)
        y = y.*(x-t(j))/(t(i+1)-t(j));
    end
end
end
