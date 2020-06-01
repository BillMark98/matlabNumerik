% mode = 1, means beta = 1/2 + sqrt(3)/6, mode = 2, means beta = 1/2 -
% sqrt(3)/6
function y = linearSDIRK(mode, A,g,n,y0,t0,tend)
if(mode == 1)
    beta = 1/2 + sqrt(3)/6;
elseif(mode == 2)
    beta = 1/2 - sqrt(3)/6;
else
    error("usage: mode 1 or 2")
end
[m,p] = size(A);
if(m ~= p)
    error("incorrect matrix size, A must be a square matrix")
end
h = (tend - t0)/n;
A_normMax = norm(A,Inf);
% if(h * beta * A_normMax > 1)
%     error("the step is too large, matrix cannot be inverted")
% end
In = eye(m);
Q = inv(In - h * beta * A);
tj = t0;
y = zeros(m,n+1);
y(:,1) = y0;
alpha = [beta, 1 - beta];
B = [beta,0;1 - 2 * beta, beta];

for j = 1 : n
    w = Q * A * y(:,j);
    k1 = w + Q * g(tj + alpha(1) * h);
    k2 = w + Q * A * h * B(2,1) * k1 + Q * g(tj + alpha(2) * h);
    y(:,j + 1) = y(:,j) + h * 1/2 * (k1 + k2);
    tj = tj + h;
end



% 
% LHS = [1 0; 0 1] - h * beta * A;
% [L,U,P] = lu(LHS)
% 
% rhs1 = A * y(:,m - 1) + g(t + alpha(1) * h);
% k1 = U\(L\(P * rhs1));
% rhs2 = A * (y(:,m-1) + h * B(2,1)*k1) + g(t + alpha(2) * h);
% k2 = U\(L\(P * rhs2));
% y(:,m) = y(:,m-1) + h * (gamma(1) * k1 + gamma(2) * k2);
