function c = Tscheby_poly(n)
if n == 0
    c = 1;
    return
end
if n == 1
    c = [1,0]';
    return
end

c0 = 1;
c1 = [1,0]';
cx = c1; 
% p_{n} = 2x * p_{n-1} - p_{n-1}
for i = 2 : n
    c0 = [0;0;c0(:,1)];
    c = 2 * conv(cx,c1) - c0;
    c0 = c1;
    c1 = c;
end
end