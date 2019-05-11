function [xk,iter] = newton( x0, tol, max_iter)
f1 = @(x,y) x^3 - 3*x*y^2 - 1;
f2 = @(x,y) 3*x^2*y - y^3;
df1x = @(x,y) 3*x^2 - 3*y^2;
df1y = @(x,y) -6*x*y;
df2x = @(x,y) 6*x*y;
df2y = @(x,y) 3*x^2 - 3*y^2;

% moegliche NS liegt in der Umgebung von drei Pkten
% (-0.5,-0.8), ( -0.5,0.8),(0.9,0.1)

x_prev = x0;
F = [f1(x_prev(1),x_prev(2));f2(x_prev(1),x_prev(2))];
A = [df1x(x_prev(1),x_prev(2)), df1y(x_prev(1),x_prev(2));
    df2x(x_prev(1),x_prev(2)),df2y(x_prev(1),x_prev(2))];

x_next = x_prev - A\F;
count = 1;
while ((norm(F) >= tol)  || (norm(x_prev - x_next) >= tol)) && (count <= max_iter)
    x_prev = x_next;
    F = [f1(x_prev(1),x_prev(2));f2(x_prev(1),x_prev(2))];
    A = [df1x(x_prev(1),x_prev(2)), df1y(x_prev(1),x_prev(2));
    df2x(x_prev(1),x_prev(2)),df2y(x_prev(1),x_prev(2))];
    x_next = x_prev - A\F;
    count = count + 1;
end

xk = x_next;
% F;
iter = count;
end
    