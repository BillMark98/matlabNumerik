function c = L2_BN(f,k)
d = zeros((k+1),1);
c = d;
int_fx = zeros((k+1),1);

for i = (k+1) : -1 : 1
    int_fx(i) = integral(@(x)(f(x) .* x.^(k+1 - i)),-1,1);
    d(i) = Legendre_poly(k+1-i)' * int_fx(i : end,1);
    d(i) = d(i) * (2*(k+1-i)+1)/2;
    c(i:end,1) = c(i:end,1) + d(i) * Legendre_poly(k+1-i);
end
end