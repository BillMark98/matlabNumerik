% Uebung 6

%% Aufgabe 2
p = @(x) x.^3 - 2*x.^2 + x;
q1 = @(x) x.^2 + 10 * x - 11;
p_q1 = @(x) x.^3 - 2*x.^2 + x - (x.^2 + 10 * x - 11);
subplot(2,1,1);
fplot(p,[-3,5],'color','blue');
hold on;
fplot(q1,[-3,5],'color','red');
subplot(2,1,2);
fplot(p_q1,[-3,5],'color','green');



%    
% u = [-1 2 3 -2 0 1 2]
% v = [2 4 -1 1]
% w = conv(u,v,'same')
% w = conv(u,v)
% 
% a = [2,1,3]'
% b = [1,0]'
% c = conv(a,b)
% 

%% test Legendre_poly and Tscheby_poly
c = Legendre_poly(4)
d = Tscheby_poly(4)
% 

%% Test L2_BN and L2c_BN
c0 = [0,1]';
c0 = [0;c0(:,1)]
f1 = @(x) (35*x.^3 - 30*x.^2 +3*x - 1)
L2_BN(f1,3)
% 
% 
f2 = @(x) (8*x.^4 - x)
L2c_BN(f2,4)


%% 
f = @(x) x.^6 - 1;
p5 = L2_BN(f,5);
p5_tild = L2c_BN(f,5);
temp = zeros(7,1);
temp(1) = 1;
% 
p5_hat = -(2.^(-5)).*Tscheby_poly(6);
p5_hat(1) = p5_hat(1) +1;
p5_hat(7) = p5_hat(7) -1;
figure(3)
x = -1:0.01:1;
plot(x,f(x),'k')
hold on
plot(x,polyval(p5,x),'r','LineWidth',3);
plot(x,polyval(p5_tild,x),'b','LineWidth',3);
plot(x,polyval(p5_hat,x),'c','LineWidth',1);
legend('f','L_2 approximation','L_{2c} approximation','L_\infty approximation')
% 
res1 = max(abs(f(x) - polyval(p5_tild,x)));
res2 = max(abs(f(x) - polyval(p5_hat,x)));
% 
norm_pmax = norm(polyval(p5_tild,x) - (x.^6 - 1),'inf')
norm_pmax = norm(polyval(p5_hat,x) - (x.^6 - 1),'inf')

%% extra tests for Aufgabe 2
x = -3 : 0.01 : 5;
p = [1,-2,1,0];
p1 = [-11/4,55/64,-1/256];
plot(x,polyval(p,x));
hold on
plot(x,polyval(p1,x));
% 
q = [1,-10,1,0];
q1 = [-43/4,55/64,-1/256];
figure
plot(x,polyval(q,x));
hold on
plot(x,polyval(q1,x));

%% plots for legendre and Tschebycheff
figure(1)
hold on
x = linspace(-1,1,101);
n = 10;
for i = 0 : n
    plot(x,polyval(Legendre_poly(i),x));
end
title('Legendre-Polynome');
hold off;

figure(2)
hold on
for i = 0 : n
    plot(x,polyval(Tscheby_poly(i),x));
end
title('Tschebyscheff-Polynome');
hold off;
