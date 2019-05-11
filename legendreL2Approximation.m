function u = legendreL2Approximation(n,f)
A = zeros((n+1),(n+1));
b = zeros((n+1),1);
for i = 1 : (n+1)
    A(i,i) = 2/(2*(i-1) + 1);
    L = @(x) ((eye(1,i) * legendre(i-1,x)));
    b_func = @(x) (L(x).* f(x));
    b(i) = integral(b_func,-1,1);
end
u = A\b;
end
