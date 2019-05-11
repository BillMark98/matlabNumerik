function [s] = Spektrum(A,tol)
n = length(A);
H = Hessenberg(A);
s = zeros(n,1);
for i = n : (-1) : 2
    [H] = QR_EW(H,tol);
    s(i) = H(i,i);
    H = H(1:(i-1),1:(i-1));
end
s(1) = H(1,1);
end
    