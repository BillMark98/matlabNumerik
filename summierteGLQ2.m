function [IQ]=summierteGLQ2(a,b,n,f)
h = (b - a)/n;
IQ = 0;
x0 = -sqrt(3/5);
x1 = 0;
x2 = sqrt(3/5);
% omega0 = 5/9;
% omega1 = 8/9;
% omega2 = 5/9;
for i = 0 : (n-1)
    x0_hat = h/2 * x0 + a + (2 * i + 1)/2 * h;
    x1_hat = h/2 * x1 + a + (2 * i + 1)/2 * h;
    x2_hat = h/2 * x2 + a + (2 * i + 1)/2 * h;
    IQ = IQ + h/18 * (5 * f(x0_hat) + 8 * f(x1_hat) + 5 * f(x2_hat));
end
end
