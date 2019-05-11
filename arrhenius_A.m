function A = arrhenius_A(K,T,E)
R = 8.314;
B = arrayfun(@(x) -1/(x*R),T);
B = E * B;
B = arrayfun(@(x) exp(x),B);
Q = B'*B; % die Matrix A'*A in der Vorlesung
A = Q\(B' * K);
end