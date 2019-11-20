function [y,err] = RK45(t0,tend,n,y0,funhandle)
h = (tend - t0)/n;
y = zeros(2,n + 1);
err = zeros(n+1,1);
err(1) = 0;
y(:,1) = y0;
t = t0;
k = zeros(2,6);
alpha = [0,1/4,3/8,12/13,1,1/2];
beta = zeros(6,6);
beta(2,1) = 1/4;
beta(3,1:2) = [3/32, 9/32];
beta(4,1:3) = [1932/2197, -7200/2197, 7296/2197];
beta(5,1:4) = [439/216, -8, 3680/513,-845/4104];
beta(6,1:5) = [-8/27,2,-3544/2565,1859/4104,-11/40];
gamma4 = [25/216,0,1408/2565,2197/4104,-1/5];
gamma5 = [16/135,0,6656/12825,28561/56430,-9/50,2/55];
for m = 2 : n+1
    for i = 1 : 6
        z = y(:, m-1);
        l = 1;
        while l < i
            z = z + h* beta(i,l) * k(:,l);
            l = l + 1;
        end
        k(:,i) = funhandle(t + alpha(i) * h, z);
    end
    y(:,m) = y(:,m-1);
    for i = 1 : 5
        y(:,m) = y(:,m) + h * gamma4(i) * k(:,i);
    end
    y5 = y(:,m-1);
    for i = 1 : 6
        y5 = y5 + h * gamma5(i) * k(:,i);
    end
    err(m) = y5(1) - y(1,m);
    t = t + h;
end
