n = 6;
f = @(t) (t  - t.^2) .* (sin(t)).^(2 * n);
F = @(x) integral(f,0,x);
x = 0:0.01:20;
y = zeros(length(x),1);
for i = 1 : length(y)
    y(i) = F(x(i));
end

plot(x,y)