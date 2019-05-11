function [B] = matrixB(n)
h = 1/n;
temp = ones(n-1,1);

A = diag(2/h^2 * temp) - diag(1/h^2 * (temp(1:(n-2))), 1) - diag(1/h^2 * (temp(1:(n-2))), -1);
R = diag(2 * temp);
R_root = sqrtm(R); % matrix square root
R_rinv = inv(R_root);
B = R_rinv * A * R_rinv;
end