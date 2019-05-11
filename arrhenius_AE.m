function [A,E] = arrhenius_AE(K,T)
m = length(T); 
R = 8.314;
B = zeros(m,2);
B(:,1) = arrayfun(@(x) 1/(-x*R),T);
B(:,2) = 1; % wollen B*[ E;log(A)] = [logK] loesen
Q = B'* B;
K = arrayfun(@(x) log(x),K);
J = Q\(B'*K);
E = J(1);
A1 = J(2);

kappa = cond(B);
cost = norm(B*J)/norm(K);
fprintf("kappa: %f\n",kappa);
fprintf("cos(theta): %f\n",cost);
A = exp(A1);
end

% [Q,R] = qr(B);
% K = arrayfun(@(x) log(x),K);
% K1 = Q'*K;
% K2 = K1(1:2,1);
% R1 = R(1:2,1:2);
% J = R1\K2;
% E = J(1);
% A1 = J(2);
% A = exp(A1);
