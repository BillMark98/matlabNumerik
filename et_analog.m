% p 194 huachengying
% result is comple i == 200
% if use log(0.73/x) will not converge 
eps = 1e-9
f_a = @(x) 26/3 * log((730)./x) - x;
f_d = @(x) -26/3 * x - 1;
x0 = 100;
xpre = x0;
for i = 1 : 30
    xnew = xpre - f_a(xpre)/(f_d(xpre));
    if(abs(xnew - xpre) < eps)
        i;
        break;
    end
    xpre = xnew;
end

xnew;

% use another method the fix point 
% x - 26/3 log ( 730/x) the unit is micro ampere
% at the fix point ~ 28.1 the derivative is 26/(3 * 28) so ok
% and for x in [10,100] . 26/3 log(730/x) is in [17.23,37.18] banach ok
% if instead use x - 26/3 log ( 0.73/x) unit is mili ampere
% the fix point is not convergent
eps = 1e-9
f_a = @(x) 26/3 * log((730)./x)
x0 = 100;
xpre = x0;
for i = 1 : 50
    xnew = f_a(xpre);
    if(abs(xnew - xpre) < eps)
        i;
        break;
    end
    xpre = xnew;
end

xnew;