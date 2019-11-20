function y = ABM3(t0,tend,n,y0,y1,y2,y3,funhandle)
h = (tend - t0)/n;
y = zeros(length(y0),n + 1);
y(:,1:4) = [y0,y1,y2,y3];
abB = [23,-16,5];
amB = [19,-5,1];
t = t0 + 4 * h;
for m = 5 : n + 1
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