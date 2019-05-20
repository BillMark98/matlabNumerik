function c = L2c_BN(f,k)
% d = zeros((k+1),1);
% c = d;
% int_fx = zeros((k+1),1);
% 
% for i = (k+1) : -1 : 1
%     int_fx(i) = integral(@(x)(f(x) .* x.^(k+1 - i) ./sqrt(1-x.^2)),-1,1);
%     d(i) = Tscheby_poly((k+1-i))' * int_fx(i : end,1);
%     if i == (k+1)
%         d(i) = d(i) / pi;
%     else 
%         d(i) = d(i)/(pi/2);
%     end
%     c(i:end,1) = c(i:end,1) + d(i) * Tscheby_poly(k+1-i);
%     
% end

w = @(x) 1./sqrt(1-x.^2); % Gewichtfunktion
c = zeros(k+1,1);
for i = 0 : k
    % the i-th tscheby polynomial
    ek = Tscheby_poly(i);
    % calculate (f,ek)/(ek,ek) 
    f_e = @(x) w(x) .* f(x) .* polyval(ek,x);
    e_e = @(x) w(x) .* polyval(ek,x).^2;
    coeff = integral(f_e,-1,1) / integral(e_e,-1,1);
    
    % convert into monom basis
    tscheb = [zeros(k - i,1);coeff * ek];
    c = c + tscheb;
     
end