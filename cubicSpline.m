function pp = cubicSpline(f,a,b)
n = length(f);
h = (b - a)/(n - 1);
m = get2ndDerivSpline(f,h);
i = 0 : (n - 1);
x = a + i * h;
L = [];
for i = 1 : (n - 1)
    ai = (m(i+1) - m(i))/(6 * h);
    bi = m(i)/2;
    ci = (f(i+1) - f(i) - m(i+1)/6 * h^2 - m(i)/3 * h^2)/h;
    di = f(i);
    L = [L;[ai bi ci di]];
end
pp = mkpp(x,L);
end
