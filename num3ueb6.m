%% b) AB4, lambda = -2
t0 = 0;
tend = 2;
y0 = 1;
lambda = -2;
f = @(t,y) lambda * y - (lambda + 1) * exp(-t);
n = 100;
y = AB4(t0,tend,n,y0,f);
timeInterval = linspace(t0,tend,n);
plot(timeInterval,y(1:length(y)-1),'r');
hold on
plot(timeInterval,exp(-timeInterval),'b');
legend('AB4','exact solution');


%% d) lambda = -2
lambda = -2;
n = 100;
h = (tend - t0)/n;
f = @(t,y) lambda * y - (lambda + 1) * exp(-t);
y = classicRK(t0,t0 + 3*h,3,y0,f);
yabm3 = ABM3(t0,tend,n,y(1),y(2),y(3),y(4),f);
timeInterval = linspace(t0,tend,n);
plot(timeInterval,yabm3(1:length(yabm3)-1),'r');
hold on
plot(timeInterval,exp(-timeInterval),'b');
legend('ABM3','exact solution');

%% d) lambda = -20
t0 = 0;
tend = 2;
y0 = 1;
lambda = -20;
n = 1000;
h = (tend - t0)/n;
f = @(t,y) lambda * y - (lambda + 1) * exp(-t);
% calculate the start value using classical Runge-Kutta method
y = classicRK(t0,t0 + 3*h,3,y0,f);
yabm3 = ABM3(t0,tend,n,y(1),y(2),y(3),y(4),f);
timeInterval = linspace(t0,tend,n);
plot(timeInterval,yabm3(1:length(yabm3)-1),'r');
hold on
plot(timeInterval,exp(-timeInterval),'b');
legend('ABM3','exact solution');

%% e) lambda = -2
% in this case, the error decreases as h decreases forr both method,
% and the error of ABM3 is smaller than that of AB4
t0 = 0;
tend = 2;
lambda = -2;
y0 = 1;
f = @(t,y) lambda * y - (lambda + 1) * exp(-t);
f_exact = @(t) exp(-t);
err = zeros(2,8);
for i = 1 : 8
    n = 2^(i + 1);
    yab4 = AB4(t0,tend,n,y0,f);
    % calculate the start value using classical Runge-Kutta method
    y = classicRK(t0,t0 + 3*h,3,y0,f);
    yabm3 = ABM3(t0,tend,n,y(1),y(2),y(3),y(4),f);
    err(1,i) = abs(yab4(n) - f_exact(2));
    err(2,i) = abs(yabm3(n) - f_exact(2));
end
semilogy(1:8,err(1,:),'r')
hold on
semilogy(1:8,err(2,:),'b')
legend('ab4','abm3')
title('two methods for ODE, when lambda = -2')
for i = 1 : 7
    order1 = log(err(1,i)/err(1,i + 1))/log(2);
    order2 = log(err(2,i)/err(2,i + 1))/log(2);
    disp1 = ['the order of convergence for AB4 when h = ' , num2str(2^(-i)), ' is ', num2str(order1)];
    disp2 = ['the order of convergence for ABM3 when h = ' , num2str(2^(-i)), ' is ', num2str(order2)];
    disp(disp1);
    disp(disp2);
end
%% e) lambda = -20
% in this case the error of both methods arise as h becomes smaller and
% then decreases as h goes further smaller.
% But in overall cases, the error of ABM3 is smaller than that of AB4
t0 = 0;
tend = 2;
lambda = -20;
y0 = 1;
f = @(t,y) lambda * y - (lambda + 1) * exp(-t);
f_exact = @(t) exp(-t);
err = zeros(2,8);
for i = 1 : 8
    n = 2^(i + 1);
    yab4 = AB4(t0,tend,n,y0,f);
    % calculate the start value using classical Runge-Kutta method
    y = classicRK(t0,t0 + 3*h,3,y0,f);
    yabm3 = ABM3(t0,tend,n,y(1),y(2),y(3),y(4),f);
    err(1,i) = abs(yab4(n) - f_exact(2));
    err(2,i) = abs(yabm3(n) - f_exact(2));
end
semilogy(1:8,err(1,:),'r')
hold on
semilogy(1:8,err(2,:),'b')
legend('ab4','abm3')
title('two methods for ODE, when lambda = -20')
for i = 1 : 7
    order1 = log(err(1,i)/err(1,i + 1))/log(2);
    order2 = log(err(2,i)/err(2,i + 1))/log(2);
    disp1 = ['the order of convergence for AB4 when h = ' , num2str(2^(-i)), ' is ', num2str(order1)];
    disp2 = ['the order of convergence for ABM3 when h = ' , num2str(2^(-i)), ' is ', num2str(order2)];
    disp(disp1);
    disp(disp2);
end