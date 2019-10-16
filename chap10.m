% A = gallery(5);
% while 1
%     clc
%     svd(A + eps*randn(5,5).*A)
%     pause(.25)
% end

% load detail
% 
% [U,S,V] = svd(X,0);
% sigma = diag(S);
% semilogy(sigma,'.')
% L = [1,20,100];
% figure(2)
% subplot(2,2,1)
% image(X)
% colormap(gray(64))
% axis image,axis off
% r = rank(X)
% title(['rank = ' int2str(r)])
% for count = 1 : 3 
%     i = L(count);
%     B = zeros(size(X));
%     for j = 1 : i
%         B = B + sigma(j) * U(:,j) * V(:,j)';
%     end
%     subplot(2,2,count + 1);
%     image(B)
%     colormap(gray(64))
%     axis image,axis off
%     title(['rank = ' int2str(i)])
% end




%% Problem

%% 10.4
n = 4;
d = ones(n-1,1);
A = diag(d,1) + diag(d,-1);
e = eig(A);
plot(-(n-1)/2 : (n-1)/2,e,'.')
figure(2)
plot(1 : n, e,'.')
% eig = 2 cos(k*pi/(n + 1)
det(A)
i = 1 : n;
ee = 2*sin(2*pi/(n+1)*i);

%% 10.5

n = 10;

figure(1)
hold on
L = linspace(0.01,0.99,100);
m = length(L);
for index = 1 : m
    t = L(index);
    A = zeros(n,n);
    for k = 1 : 10
        for j = 1 : 10
            A(k,j) = 1/(k - j + t);
        end
    end
    e = eig(A);
    plot(real(e),imag(e), '.')
    xlim([-10,10]);
    ylim([-10,10]);
end

%% 10.6
A = gallery(5)
[X,w] = eig(A)
det(X)
v = condeig(A)

%% 10,7
R = sym(rosser)
e = eig(R)
% e is symbols, so have to convert to double , so that we can sort
[ignore,k] = sort(double(e))
e = e(k)
% poly removed, use charpoly
p = charpoly(R)
pp = poly2sym(p)
syms x
f = factor(pp,x)
pretty(f)

e = eig(sym(rosser))
r = eig(rosser)
double(e)-r
double(e-r)

%%
n = 10
A = diag(ones(n-1,1),-1) + diag(1,n-1)

eigsvdgui(A,'eig')
