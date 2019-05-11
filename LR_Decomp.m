function [L,U,p] = LR_Decomp(A)
% [m,n] = size(A);
% if m ~= n
%     disp('size mismatch');% not a square matrix, LR-decomposition not exist
%     return
% end
[m,m] = size(A);
LR = A;
p = 1:m;
for j = 1:m
    max = j;
    for k = j:m % traverse the row for the j-th column to find the
        % maximum
        if abs(LR(k,j)) > abs(LR(max,j))
            max = k;
        end
    end
    if abs(LR(max,j)) < 1e-9
        fprintf("the value %f\n",LR(max,j))  % det(A) = 0
        disp('Matrix singular')
        return
    end
    for k = 1:m  % calculate the permutation 
        % to get the permutation matrix P use P = eye(m); P = P(:,p)
        if p(k) == j
            p(k) = max;
        elseif p(k) == max
            p(k) = j;
        end
    end
    permutemp = [1:m];
    permutemp(j) = max;
    permutemp(max) = j;
    LR = LR(permutemp, :); %permute the matrix
    for i = (j+1):m
        LR(i,j) = LR(i,j)/LR(j,j); %calculate the corresponding 
                                               %column of the L matrix
        for t = (j+1) : m
            LR(i,t) = LR(i,t) - LR(j,t)*LR(i,j); % vector version : LR(i,(j+1):end) = LR(i,(j+1):end)-LR(i,j)*LR(j,(j+1):end)
        end
    end
end
L = tril(LR,-1) + eye(m,m);
U = triu(LR);
end