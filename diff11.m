syms K x
% S = vpasolve([K/(1 + (K/41 - 1)* exp(-x * 30)) == 56, K/(1 + (K/41 - 1)* exp(-x * 100)) == 78],[K,x])
% S = vpasolve([log(56/41) - log((K - 56000000)/(K - 41000000)) == x * 30, log(78/41) - log((K - 78000000)/(K - 41000000)) == x * 100],[K,x])
S = vpasolve([log(56/41) - log((K - 56)/(K - 41)) == x * 3, log(78/41) - log((K - 78)/(K - 41)) == x * 10],[K,x])
fprintf("S.K is %e, S.x is %e",S.K,S.x)

k = S.K
alpha = S.x
% res1 = k/(1 + (k/41 - 1)* exp(-alpha * 3))
% res2 = k/(1 + (k/41 - 1)* exp(-alpha * 10))


res30 = log(56/41) - log((k - 56)/(k - 41))
res31 = alpha * 3

res41 = log(78/41) - log((k - 78)/(k - 41))
res42 = alpha * 10