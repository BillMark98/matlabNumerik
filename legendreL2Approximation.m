function u = legendreL2Approximation(n,f)
A = zeros((n+1),(n+1));
b = zeros((n+1),1);
for i = 1 : (n+1)
    A(i,i) = 2/(2*(i-1) + 1);
    % legendre(N,x) returns a (N + 1) * length(x) matrix ,denote by A
    % each line is the m-th derivative of the n-th Legendre polynomial 
    % evaluated at point x(k) k = 1..end
    % A_{1,-} = Pn(x(1)), Pn(x(2)),...,Pn(x(end))
    % A_{2,-} = Pn'(x(1)), Pn'(x(2)), ... , Pn'(x(end))
    % A_{3,-} = Pn''(x(1)), ......
    L = @(x) ((eye(1,i) * legendre(i-1,x)));
    b_func = @(x) (L(x).* f(x));
    b(i) = integral(b_func,-1,1);
    u(i) = b(i)/A(i,i);
end
fprintf("condition of matrix A legendreapprox with n = %d is %f\n",n,cond(A));
end
