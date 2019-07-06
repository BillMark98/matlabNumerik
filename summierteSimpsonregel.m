function [IS] = summierteSimpsonregel(a,b,n,f)
h = (b - a)/n;
IS = 0;
for i = 0 : (n-1)
    xi = a + i*h;
    IS = IS + 1/6 * h * (f(xi) + 4 * f(xi + 1/2 * h) + f(xi + h));
end
end