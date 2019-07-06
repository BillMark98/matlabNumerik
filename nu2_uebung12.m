% f = @(x) exp(x)./(x.^2 + 0.5)
% sum = 1/2 * 1/6*(f(0) + 4*f(1/4) + 2*f(1/2)+ 4 * f(3/4) + 2*f(1) + 4 * f(5/4) + 2 * f(3/2) + 4 * f(7/4) + f(2))
% integral(f,0,2)
% 
% f1 = @(x) 1./(1 + x.^2)
% h = 1/2;
% sum1 = h/2*(f1(h*(1/2-sqrt(3)/6)) + f1(h * (1/2 + sqrt(3)/6)) + f1( 1/2 + h*(1/2-sqrt(3)/6)) +...
%     f1(1/2 + h * (1/2 + sqrt(3)/6)))
% Ef = h^4/4320 * 24
% integral(f1,0,1)
a = 0;b = 1;
N = 10;
F = @(x) 1./(1 + x)
err1 = zeros(N+1,1);
err2 = zeros(N+1,1);
index = 0 : N;
nn = 2.^index;
h = 1./nn;
for i = 0 : N
    n = 2^i;
    IS = summierteSimpsonregel(a,b,n,F);
    IQ = summierteGLQ2(a,b,n,F);
    err1(i+1) = abs(IS - log(2))/log(2);
    err2(i+1) = abs(IQ - log(2))/log(2);
end
loglog(h,err1,'-s');
grid on;
hold on;
loglog(h,err2,'-s');
legend('simpson','gaus')
