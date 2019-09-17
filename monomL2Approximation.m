function [u] = monomL2Approximation(n,f)
A = zeros((n+1),(n+1));
b = zeros((n+1),1);
for i = 1 : (n+1)
%     for j = i : (n+1)
%         A(i,j) = integral(@(x) x.^(i-1).*x.^(j-1),-1,1); % should use .^ .* or else would 
%         % throw error message use .^ .* for elementwise multi
%         if i ~= j
%             A(j,i) = A(i,j);
%         end
%     end
    b(i) = integral(@(x) x.^(i-1) .* f(x), -1,1);
    
    % another version for calculation of A
    % since index begins with 1, which means, e.g A(1,1) = (p_0 , p_0) = 2
    % = 2/(2 * 1 - 1)  A(i,i) = (x^(i-1) , x^(i-1)) = 1/(2i-1) * 2
    A(i,i) = 2/(2*i - 1);
    for j = (i+1):(n+1)
        if (mod(i+j,2) == 0)
             A(i,j) = 2/(i+j-1);
             A(j,i) = 2/(i+j-1);
        end
    end
    
end
fprintf("condition of matrix A monomapprox with n = %d is %f\n",n,cond(A));
u = A\b;
end