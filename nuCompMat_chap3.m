% 18
t = 1900 : 10 : 2000;
V = vander(t);
co1 = cond(V)
mu = mean(t);
sigma = std(t);
s = (t - mu)/sigma;
V = vander(s);
co2 = cond(V)