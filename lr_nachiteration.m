function [x,r] = lr_nachiteration(A,LR,p,b,k)
L = tril(LR,-1);
[m,n] = size(L);
if m ~= n
    disp('size mismatch')
    return
end
L = L + eye(m); % calculate L
R = triu(LR); 
P = eye(m);
P = P(p,:); % Calculate the permutation matrix
x = zeros(m,1);
x_prev = x;
y = ones(m,1);

b = P*b;
r = b;
% Calculate the y through the equation Ly = Pb
for count = 1 : k
    for i = 1 : m
        sum = 0;
        for j = 1 : (i-1)
            sum = sum + y(j) * L(i,j);
        end
        y(i) = (r(i) - sum);
    end
    % Calculate x through the equation Rx = y
    for i = m : (-1) : 1
        sum = 0;
        for j = m : (-1) : (i + 1)
            sum = sum + x_prev(j) * R(i,j);
        end
        x_prev(i) = (y(i) - sum) / R(i,i);
    end
   
    x = x + x_prev;
    r = b - A*x;
end
end
    