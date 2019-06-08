x = 1;
f = @(x) x.^(1./x);
Df = @(x) x.^(1./x) .* (1 - log(x))./(x.^2);
DDf=@(x) Df(x).*((-3./x.^3 + 2 * log(x)./(x.^3)) + ((1 - log(x))./(x.^2)).^2);
DDDf = @(x) Df(x).*(((1 - log(x))./(x.^2)).^3 + 3 * ((1 - log(x))./(x.^2)) .* (-3./x.^3 + 2 * log(x)./(x.^3))...
    + (-6*log(x)./(x.^4) + 11/x.^4));

i = 1:10;
h = 10.^(-i);
errp = abs(vorwaertsDifferenz(h,x)-Df(x));
errc = abs(zentraleDifferenz(h,x)-Df(x));
subplot(2,1,1)
loglog(h,errp,'-*');
% optimal h at 10^(-8)
% slope 1
% calculate the theoretic value
% h = sqrt(4 * eps * || f ||_infty / || f''||_infty) here since h is small
% || f ||_infty is f(1) and || f'' ||_infty is |f''(1)|
hp_theory = sqrt(4*eps*abs(f(1))/abs(DDf(1)));
legend('errp');
subplot(2,1,2)
loglog(h,errc,'+');
legend('errc');
%optimal h at 10^(-7)
%slope 2

% calculate the theoretic value
% h = (3 * eps * || f ||_infty / || f'''||_infty)^(1/3) here since h is small
% || f ||_infty is f(1) and || f'' ||_infty is |f''(1)|
hc_theory = (3 * eps * abs(f(1))/abs(DDDf(1)))^(1/3);