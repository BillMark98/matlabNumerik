function c = L2_BN(f,k)
% d = zeros((k+1),1);
% c = d;
% int_fx = zeros((k+1),1);
% 
% for i = (k+1) : -1 : 1
%     int_fx(i) = integral(@(x)(f(x) .* x.^(k+1 - i)),-1,1);
%     d(i) = Legendre_poly(k+1-i)' * int_fx(i : end,1);
%     d(i) = d(i) * (2*(k+1-i)+1)/2;
%     c(i:end,1) = c(i:end,1) + d(i) * Legendre_poly(k+1-i);
% end

c = zeros(k+1,1);
for i = 0 : k
    legend_i = Legendre_poly(i);
    f_k = @(x) f(x).* polyval(legend_i,x);
    e_k = @(x) polyval(legend_i,x).^2;
    
    % the coeff of the i-th legendre polynomial
    coeff = integral(f_k,-1,1)/integral(e_k,-1,1);
    
    % change into the monomial base
    legend_cof = [zeros(k-i,1); coeff * legend_i];
    c = c + legend_cof;
end