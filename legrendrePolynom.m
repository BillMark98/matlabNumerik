% given u, and x calculate the polynomial evaluted at x(k) k = 1..end
% where the polynomial is given by sum(u(i) * P_{i-1}(x), i = 1..length(u))
% where P_{i-1} is the i-th Legendre polynomial ( degree = i-1)
function y = legrendrePolynom(u,x)
n = length(u);
m = length(x);
y = zeros(m,1);
% u = col2RowTransform(u);
% for i = 1 : m
%     sum = 0;
%     for j = 0 : (n-1)
%         sum = sum + eye(1,(j+1)) * legendre(j,x(i)) * u(j+1);
%     end
%     y(i) = sum;
% end
sum = zeros(1,m);
for j = 0 : (n - 1)
    sum = sum + eye(1,(j + 1)) * legendre(j,x) * u(j + 1);
end
y = sum';
end

function u = col2RowTransform(x)
[m,n] = size(x);
if (m > n)
    u = x';
else
    u = x;
end
end
