function [u] = monomL2Approximation(n,f)
A = zeros((n+1),(n+1));
b = zeros((n+1),1);
for i = 1 : (n+1)
    for j = i : (n+1)
        A(i,j) = integral(@(x) x.^(i-1).*x.^(j-1),-1,1); % should use .^ .* or else would 
        % throw error message use .^ .* for elementwise multi
        if i ~= j
            A(j,i) = A(i,j);
        end
    end
    b(i) = integral(@(x) x.^(i-1) .* f(x), -1,1);
end
u = A\b;
end