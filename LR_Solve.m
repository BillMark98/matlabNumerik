function x = LR_Solve(LR,p,b)
L = tril(LR,-1);
[m,n] = size(L);
if m ~= n
    disp('size mismatch')
    return
end
L = L + eye(m); % Calculate y
R = triu(LR); % Calculate R
P = eye(m);
P = P(p,:); % Calculate the permutation matrix
x = ones(m,1);
y = ones(m,1);
b = P*b;
% Calculate the y through the equation Ly = Pb
for i = 1 : m
    sum = 0;
    for j = 1 : (i-1)
        sum = sum + y(j) * L(i,j);
    end
    y(i) = (b(i) - sum);
end
% Calculate x through the equation Rx = y
for i = m : (-1) : 1
    sum = 0;
    for j = m : (-1) : (i + 1)
        sum = sum + x(j) * R(i,j);
    end
    x(i) = (y(i) - sum) / R(i,i);
end
end
