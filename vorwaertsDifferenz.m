function y = vorwaertsDifferenz(h,x)
f = @(x) x.^(1./x);
y = (f(x+h) - f(x))./h;
end