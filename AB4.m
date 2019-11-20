function y = AB4(t0,tend,n,y_init,funhandle)
h = (tend - t0)/n;
y0 = y_init(:,1);
y = zeros(length(y0),n + 1);
% initialize the y
y(:,1:4) = y_init(:,1:4);
% the coeffcients of bi of AB method, the factor 1/24 will be multiplied
% after summing up all the values, see the code below.
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
