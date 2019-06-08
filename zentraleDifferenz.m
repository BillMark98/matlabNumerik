% calculate for f = 1/(1+x^2) the central difference
% D_h^{c}(x) = (f(x+h)-f(x-h))/2h
function y = zentraleDifferenz(h,x)
f = @(x) x.^(1./x);
y = (f(x+h) - f(x-h))./(2*h);
end