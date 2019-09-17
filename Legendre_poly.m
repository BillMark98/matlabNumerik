% given n, give back a vector c
% such that the n-th Legendre Polynomial (degree n) is given by
% c(1)x^n + c(2)x^(n-1).... + c(n + 1)
function c = Legendre_poly(n)
if n == 0
    c = 1;
    return
end
if n == 1
    c = [1,0]';
    return
end

c0 = 1;
c1 = [1,0]';
cx = c1; 
% p_{n} = (2*n-1)/(n) x * p_{n-1} - (n-1)/(n)*p_{n-2}
for i = 2 : n
    c0 = [0;0;c0(:,1)];
    % conv(cx, c1) gives back the coeff of the polynomial of x * p_{n-1}
    % e.g [1,0,1] --- x^2 + 1, [1,0] ---- x
    % conv([1,0,1],[1,0]) = [1,0,1,0] --- x^3 + x
    c = (2*i-1)/i * conv(cx,c1) - (i-1)/i * c0;
    c0 = c1;
    c1 = c;
end
end