function [H] = Hessenberg(A)
n = length(A);
if n <= 2
    H = A
    return
end
for i = 1:n-2 
    v = A((i+1):n,i);
    v(1) = v(1) + sign(v(1))*sqrt(v'*v);
    Qvh = eye(n-i) - 2 * v * v'/(v'*v);
    % construct the block diagonal matrix
    % i.e Qv = diag(eye(i),Qvh)
    Qv = blkdiag(eye(i),Qvh);
    A = Qv*A*Qv;
end
H = A;
end