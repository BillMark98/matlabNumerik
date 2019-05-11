function y = legrendrePolynom(u,x)
n = length(u);
m = length(x);
y = zeros(m,1);
for i = 1 : m
    sum = 0;
    for j = 0 : (n-1)
        sum = sum + eye(1,(j+1)) * legendre(j,x(i)) * u(j+1);
    end
    y(i) = sum;
end