function y = AB4(t0,tend,n,y0,funhandle)
h = (tend - t0)/n;
y_init = classicRK(t0,t0 + 3*h,3,y0,funhandle);
y = zeros(length(y0),n + 1);
% err = zeros(n+1,1);
y(:,1) = y0;
y(:,2:4) = y_init(:,2:4);
bl = [-9, 37, -59,55];
t = t0;
for m = 5 : n + 1
    z = zeros(length(y0),1);
    for i = 1 : 4
        z = z + bl(i) * funhandle(t + (i-1) * h, y(:,m - (5 - i)));
    end
    y(:,m) = y(:,m-1) + h/24 * z;
    t = t + h;
end
