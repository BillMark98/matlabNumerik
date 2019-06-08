x = 1;
Df = @(x) x.^(1./x) .* (1 - log(x))./(x.^2);
i = 1:10;
h = 10.^(-i);
errp = abs(vorwaertsDifferenz(h,x)-Df(x));
errc = abs(zentraleDifferenz(h,x)-Df(x));
subplot(2,1,1)
loglog(h,errp,'-*');
% optimal h at 10^(-8)
% slope 1
legend('errp');
subplot(2,1,2)
loglog(h,errc,'+');
legend('errc');
%optimal h at 10^(-7)
%slope 2