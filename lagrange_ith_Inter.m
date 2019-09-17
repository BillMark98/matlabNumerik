% calculate the i'th lagrange polynomial of degree n for the points in t
% at point t
function prod = lagrange_ith_Inter(i,n,ti,t)
    prod = 1;
    for index = 1 : (n+1)
        if index ~= (i + 1)
            prod = prod * (t - ti(index))/(ti(i + 1) - ti(index));
        end
    end
end