% given u, and x calculate the polynomial evaluted at x(k) k = 1..end
% where the polynomial is given by sum(u(i) * x^(i-1), i = 1..length(u))
function y = monomPolynom(u,x)
n = length(u);
m = length(x);
y = zeros(m,1);
for i = 1 : m
    sum = 0;
    for j = 0 : (n-1) % calculate u(1) * x(i)^0 + u(2) * x(i) + u(3) * x(i)^2...
        sum = sum + u(j+1) * (x(i))^j;
    end
    y(i) = sum;
end
end

    