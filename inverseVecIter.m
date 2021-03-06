function [lambda, x,k,v_new] = inverseVecIter(A,x_j,lambda_j,tol,lambda0,k0,v)
x_j = 1/sqrt(x_j'*x_j) * x_j;
n = size(A);
D = eye(n);
z = (A - lambda_j * D)\x_j;
lambda_new = lambda_j + 1/(z' * x_j);
x_j = 1/sqrt(z'*z) * z;
k = k0 + 1;
v_new = v;
v_new(k) = log(abs(lambda_new - lambda0));
while(abs(lambda_new - lambda_j) >= tol)
    lambda_j = lambda_new;
    z = (A - lambda_j * D)\x_j;
    lambda_new = lambda_j + 1/(z' * x_j);
    x_j = 1/sqrt(z'*z) * z;
    k = k + 1;
    v_new(k) = log(abs(lambda_new - lambda0));
end
lambda = lambda_new;
x = z;