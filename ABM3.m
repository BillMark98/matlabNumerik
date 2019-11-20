function y = ABM3(t0,tend,n,y_init,funhandle)
h = (tend - t0)/n;
y0 = y_init(:,1);
y = zeros(length(y0),n + 1);
% the yi stands for the initial value at time t0 + jh;
y(:,1:3) = y_init(:,1:3);
% coeff of AB method with k = 3, the factor 12 will be multiplied after
% summing up all the values, see the code below
abB = [23,-16,5];
% the f(t_{j + 1},y_{j + 1}) is calucalted before the loop, so here only
% three coeff for the AM method k = 4;
amB = [19,-5,1];
% time begins at t0 + 2h, corresponding to the first z calculated,
% funhandle(t,y_j)
t = t0 + 2 * h;
for m = 4 : n + 1
    yj0 = y(:,m - 1);
    z = zeros(length(y0),1);
    for i = 1 : 3
        z = z + abB(i) * funhandle(t - (i-1) * h, y(:,m - i));
    end
    yj0 = yj0 + h/12 * z;
    z = 9 * funhandle(t + h, yj0);
    for i = 1 : 3
        z = z + amB(i) * funhandle(t - (i - 1) * h, y(:,m-i));
    end
    y(:,m) = y(:,m-1) + h/24 * z;
    t = t + h;
end