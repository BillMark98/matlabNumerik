function [x,res] = houseQR_Solve(QR,b)

[m,n] = size(QR);
if m < n
    disp('m < n keine Ausgleichsproblem');
    return
end
% index_end = min(m,n);
for i = 1 : n
    v = QR(i:m,i);
    v(1) = 1;
    beta = v'*v;
    z = b(i:m);
    b(i:m) = z - 2*(v'*z)/beta*v;
end
x = zeros(n,1);
for j = n : -1 : 1
    sum = 0;
    for k = n : -1 : j
        sum = sum + QR(j,k)*x(k);
    end
    x(j) = (b(j) - sum) / QR(j,j);
end
res = norm(b((n+1):end,1));
end