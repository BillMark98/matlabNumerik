format long
tau1 = 6.6e-2;
tau2 = 7.5e4;
tau3 = 1.6e3;

lambda1 = log(2)/tau1;
lambda2 = log(2)/tau2;
lambda3 = log(2)/tau3;

lambda_min = min([lambda1,lambda2,lambda3]);
lambda_max = max([lambda1,lambda2,lambda3]);
f = @(y) [-lambda1 * (y(1)); lambda1 * y(1) - lambda2 * y(2);
    lambda2 * y(2) - lambda3 * y(3)];

t0 = 0;
tend = 100;

y0 = [1;0;0];
hmin = min( [1/lambda2, 1/lambda3 , 2/lambda1])
% in this case since the first component starts from 1
% by induction y1_k = (1 - h * lambda1)^k (the k-th calculated first component
% of y) so the step size should be lower
% than min( [1/lambda2, 1/lambda3 , 2/lambda1])
% for the third component in the formula
% | 1 - h * lambda1 | <= 1 -->  -1 <= h * lambda1 - 1<= 1 
%  --> 0 <=  h * lambda1 <= 2

% but actually, we know from the lecture, that the stability region of
% euler is (-2,0), so |h * lambda| <= 2 ....


% for i = 9: 11
%     figure
%     n = 2^i;
%     h = (tend - t0)/n;
%     y = explicitEuler(f,y0,t0,tend,n);
%     hold on;
% %     plot(linspace(t0,tend,n + 1),y(1,:))
% %     plot(linspace(t0,tend,n+1),y(2,:))
% %     plot(linspace(t0,tend,n + 1), y(3,:))
% %     legend("y1","y2","y3")
% %     title(['step size:',num2str(h)])
%     
%     
%     
%     subplot(3,1,1)
%     plot(linspace(t0,tend,n + 1),y(1,:))
%     legend("y1")
%     subplot(3,1,2)
%     plot(linspace(t0,tend,n+1),y(2,:))
%     legend("y2")
%     subplot(3,1,3)
%     plot(linspace(t0,tend,n + 1), y(3,:))
%     legend("y3")
%     sgtitle(['step size:',num2str(h)])
% end

% we see that if h > 0.1904, then the solution explodes
H = [0.2, 0.1905,0.1904,0.1903, 0.1901,0.18];
for i = 1 : length(H)
    h = H(i);
    figure
    hold on;
    n = ceil((tend - t0)/h);
    y = explicitEuler(f,y0,t0,tend,n);
    subplot(3,1,1)
    plot(linspace(t0,tend,n + 1),y(1,:))
    legend("y1")
    subplot(3,1,2)
    plot(linspace(t0,tend,n+1),y(2,:))
    legend("y2")
    subplot(3,1,3)
    plot(linspace(t0,tend,n + 1), y(3,:))
    legend("y3")
    sgtitle(['step size:',num2str(h)])
end


%% b

format long
tau1 = 6.6e-2;
tau2 = 7.5e4;
tau3 = 1.6e3;

lambda1 = log(2)/tau1;
lambda2 = log(2)/tau2;
lambda3 = log(2)/tau3;

lambda_min = min([lambda1,lambda2,lambda3]);
lambda_max = max([lambda1,lambda2,lambda3]);

t0 = 0;
tend = 100;
y0 = [1;0;0];

% since in this case we can write down the formula for the implicit euler
% by basic algebra we have
% y1_{j + 1} = 1/(1 + h * lambda1) * y1_{j}
% y2_{j + 1} = 1/(1 + h * lambda2) * (y2_{j} + h * lambda1 * 1/(1 + h *
% lambda1) * y1_{j})
% y3_{j + 1} = 1/(1 + h*lambda3) * (y3_j + h * lambda2 * y2_{j+1})
% we know that h can be any positive integer
H = [0.5, 0.2,0.1904];
for index = 1 : length(H)
    h = H(index);
    h = H(index);
    figure
    hold on;
    n = ceil((tend - t0)/h);
    y = zeros(3,n+1);
    y(:,1) = y0;
    for i = 2 : (n + 1)
        y(1,i) = 1/(1 + h * lambda1) * y(1,i-1);
        y(2,i) = 1/( 1 + h * lambda2) * (y(2,i - 1) + h * lambda1 * y(1,i));
        y(3,i) = 1/(1 + h * lambda3) * (y(3,i-1) + h * lambda2 * y(2,i));
    end
    subplot(3,1,1)
    plot(linspace(t0,tend,n + 1),y(1,:))
    legend("y1")
    subplot(3,1,2)
    plot(linspace(t0,tend,n+1),y(2,:))
    legend("y2")
    subplot(3,1,3)
    plot(linspace(t0,tend,n + 1), y(3,:))
    legend("y3")
    sgtitle(['step size:',num2str(h)])
end






%% A26
alpha = @(x) (5/24 - 5/8 * exp(1i * x));
beta = @(x) (1/12 - 13/12 * exp(1i * x));
gamma = @(x) (exp(2i * x) - exp(1i * x));

f1 = @(x) (- beta(x)./(2*alpha(x)) + sqrt(beta(x).^2/(4 * alpha(x).^2) - gamma(x)./alpha(x)));
f2 = @(x) (- beta(x)./(2*alpha(x)) - sqrt(beta(x).^2/(4 * alpha(x).^2) - gamma(x)./alpha(x)));

x = linspace(-pi, pi, 100);
plot(real(f1(x)),imag(f1(x)),'bo')
hold on

plot(real(f2(x)),imag(f2(x)),'r*')
