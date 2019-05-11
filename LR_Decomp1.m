function [result,L,R,P,p] = LR_Decomp1(A)
[m,n] = size(A);
if m ~= n
    disp('size mismatch');% not a square matrix, LR-decomposition not exist
    return
end
result = A;
p = [1:m];
for j = 1:m
    max = j;
    for k = j:m % traverse the row for the j-th column to find the
        % maximum
        if abs(result(k,j)) > abs(result(max,j))
            max = k;
        end
    end
    if abs(result(max,j)) < 1e-9
        fprintf("the value %f\n",result(max,j))  % det(A) = 0
        disp('Matrix singular')
        return
    end
    for k = 1:m  % calculate the permutation matrix
        if p(k) == j
            p(k) = max;
        elseif p(k) == max
            p(k) = j;
        end
    end
    permutemp = [1:m];
    permutemp(j) = max;
    permutemp(max) = j;
    result = result(permutemp, :); %permute the matrix
    for i = (j+1):m
        result(i,j) = result(i,j)/result(j,j); %calculate the corresponding 
                                               %column of the L matrix
        for t = (j+1) : m
            result(i,t) = result(i,t) - result(j,t)*result(i,j);
        end
    end
end
L = result;
R = result;
for i = 1 : m
    L(i,i) = 1;
    for j = (i+1) : m
        L(i,j) = 0;
    end
    for j = 1 : (i - 1)
        R(i,j) = 0;
    end
    
end
P = eye(m);
P = P(p,:);
                    
                        
                        
                
end