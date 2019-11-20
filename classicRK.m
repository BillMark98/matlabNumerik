function y = classicRK(t0,tend,n,y0,funhandle)
h = (tend - t0)/n;
y = zeros(length(y0),n + 1);
y(:,1) = y0;
k = zeros(length(y0),4);
% the coeff of the alphas
alpha = [0,1/2,1/2,1];
% the coeff of the gammas
gamma = [1/6, 1/3,1/3,1/6];
t = t0;
for m = 2 : n + 1
    k(:,1) = funhandle(t,y(:,m-1));
    for i = 2 : 4
        z = k(:,(i-1));
        k(:,i) = funhandle(t + alpha(i) * h, y(:,m-1) + h * alpha(i) * z);
    end
    y(:,m) = y(:,m-1);
    for i = 1 : 4
        y(:,m) = y(:,m) + h * gamma(i) * k(:,i);
    end
    t = t + h;
end
            