function [Ak] = QR_EW(H,tol)
n = length(H);
count = 1;
while ((abs(H(n,n-1)) >= tol) || (count < 1000))
    m = H(n,n);
    [Q,R] = qr(H - m *eye(n));
    H = R*Q + m*eye(n);
    count = count + 1;
end
Ak = H;
end