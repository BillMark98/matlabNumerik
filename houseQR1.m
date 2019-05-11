function res = houseQR1(A)
[m,n] = size(A);
res = A;
for i = 1 : n
    e = zeros(m-i+1,1); % the subvector which spans the line on which we want to project
                        % our column vector of A.
    e(1) = 1;
    scale = sqrt(res(i:m,i)'*res(i:m,i)); % the length of the vector
    if (sign(res(i,i)) == 0)
        v = res(i:m,i) - scale * e; % res(i,i) = 0 , we choose here -1 instead
        pivot = scale;
    else
        v = res(i:m,i) + sign(res(i,i))* scale * e;
        pivot = scale * -1 * sign(res(i,i)); % the new res(i,i) value
    end
    v = v/v(1); % normalize the vetor v s.t. v(1) = 1
    v_len = v'*v;
    
    res(i:m,i) = v;
    res(i,i) = pivot;
    % calculate the product:  Qv * A
    for j = (i+1):n
        vec_old = res(i:m,j);
        res(i:m,j) = vec_old - 2/v_len *(v'*vec_old) * v;
    end
end
end