%% Aufgabe 28
n = 40;
h = 1/n;
v1 = [2*ones(n/2-2,1);12/11;21/110;1/5 * ones(n/2-2,1)];
v2 = [-1 * ones(n/2-2,1);-1/11;-0.1 * ones(n/2-2,1)];
A = diag(v1) + diag(v2,-1) + diag(v2,1)
b = 0.1 * ones(n-2,1)
G = 1/(h^2) * A;
u = G\b;
% calculate u_{n/2}
u_middle = 1/11 * u(n/2) + 10/11 * u(n/2 - 1);
u = [u(1:n/2 - 1,1);u_middle;u(n/2:end,1)];
% add the two values at the side 
u = [0;u;0];
x = linspace(0,1,n + 1);
plot(u,x)

% use the upwind difference to approximate first order derivative
q = [0;diff(u)];
q = [10 * q(1:n/2+1,1);q(n/2+2:end,1)];

plot(x,u)
hold on
plot(x,q)
legend('u(x)','q(x)')

%% another version
figure
n = 40;
h = 1/n;
v1 = [20*ones(n/2-1,1);11;2 * ones(n/2-1,1)];
v2 = [-10 * ones(n/2-1,1);-1 * ones(n/2-1,1)];
A = diag(v1) + diag(v2,-1) + diag(v2,1);
b = 1 * ones(n-1,1);
b(n/2 + 1,1) = 0;
G = 1/(h^2) * A;
u = G\b;
% add the two values at the side 
u = [0;u;0];
x = linspace(0,1,n + 1);
plot(u,x)

% use the upwind difference to approximate first order derivative
q = [0;diff(u)];
q = [10 * q(1:n/2+1,1);q(n/2+2:end,1)];

plot(x,u)
hold on
plot(x,q)
legend('u(x)','q(x)')

% Aufgabe 29
%% n = 10
n = 10;
h = 1/(n + 1);
a = 40;
f = @(x) (exp(a * x) - 1)./(exp(a) - 1);
A = 2 * eye(n) - diag(ones(n-1,1),-1) - diag(ones(n-1,1),1);
B = diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
b = zeros(n,1);
b(n) = a/(2 * h) - 1/(h^2);
% the final matrix
G = -1/(h^2) * A - a/(2 * h) * B;
% the solution y
y = G\b;
% add the two values at the side 
y = [0;y;1];
x = linspace(0,1,n+2);
plot(x,y,'o')
hold on
plot(x,f(x))
legend('discretized solution','exact solution')

%% n = 25
n = 25;
h = 1/(n + 1);
a = 40;
f = @(x) (exp(a * x) - 1)./(exp(a) - 1);
A = 2 * eye(n) - diag(ones(n-1,1),-1) - diag(ones(n-1,1),1);
B = diag(ones(n-1,1),1) - diag(ones(n-1,1),-1);
b = zeros(n,1);
b(n) = a/(2 * h) - 1/(h^2);
G = -1/(h^2) * A - a/(2 * h) * B;
y = G\b;
% add the two values at the side 
y = [0;y;1];
x = linspace(0,1,n+2);
plot(x,y,'o')
hold on
plot(x,f(x))
legend('discretized solution','exact solution')

%% c) n = 10
n = 10;
h = 1/(n + 1);
a = 40;
f = @(x) (exp(a * x) - 1)./(exp(a) - 1);
A = 2 * eye(n) - diag(ones(n-1,1),-1) - diag(ones(n-1,1),1);
B = eye(n) - diag(ones(n-1,1),-1);
b = zeros(n,1);
b(n) = - 1/(h^2);
G = -1/(h^2) * A - a/(h) * B;
y = G\b;
% add the two values at the side 
y = [0;y;1];
x = linspace(0,1,n+2);
plot(x,y,'o')
hold on
plot(x,f(x))
legend('discretized solution','exact solution')

%% c) n = 25
n = 25;
h = 1/(n + 1);
a = 40;
f = @(x) (exp(a * x) - 1)./(exp(a) - 1);
A = 2 * eye(n) - diag(ones(n-1,1),-1) - diag(ones(n-1,1),1);
B = eye(n) - diag(ones(n-1,1),-1);
b = zeros(n,1);
b(n) = - 1/(h^2);
G = -1/(h^2) * A - a/(h) * B;
y = G\b;
% add the two values at the side 
y = [0;y;1];
x = linspace(0,1,n+2);
plot(x,y,'o')
hold on
plot(x,f(x))
legend('discretized solution','exact solution')