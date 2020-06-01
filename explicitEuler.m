function y = explicitEuler(f,y0,t0,tend,n)
m = length(y0);
h = (tend - t0)/n;
y = zeros(m,n+1);
y(:,1) = y0;
for i = 2 : (n + 1)
    y(:,i) = y0 + f(y0) * h;
    y0 = y(:,i);
end
end


function [y] = ImplEuler(t0,te,n,y0,A)
h = (te - t0)/n;
y = zeros(3,n+1);
y(:,1) = y0;
Atilde = I - h * A;
for i = 2 : n + 1
    y(:,i) = linsolve(Atilde,y(:,i-1));
end
end