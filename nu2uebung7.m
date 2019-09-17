%% A1
f = @(x) 1./(1 + x);
% p1 = [-1/2, 1/4 + sqrt(2)/2];
q = @(x) f(x) - (-1/2 *x + 1/4 + sqrt(2)/2);
fplot(q,[0,1]);
v = [q(0),q(sqrt(2)-1),q(1)]

%% A3 d
L = [0,0.2,-1 + sqrt(2)];

prod = zeros(3,10);
for index = 1 : length(L)
    t = L(index);
    count = 1;
    for n = 1 : 50 : 451;
        i = 0 : 1 : n;
        ti = 2 * i/n - 1;
        prod(index,count) = lagrangeInter(n,ti,t);
        count = count + 1;
    end
    
end  
% % 
% % 
% 
% % another method
% 
% L = [0,0.2,-1 + sqrt(2)];
% 
% prod = zeros(3,10);
% n = 1 : 50 : 451;
% tests = max(size(n));
%     for count = 1 : tests
%         i = 0 : 1 : n(count);
%         ti = 2 * i/n(count) - 1;
%         for deg = 0 : n(count)
%             prod(:,count) = prod(:,count) + abs(lagrangePolynomial(deg,ti,L'));
%         end
%     end
% 
% 

% prod = log(prod);
% subplot(1,3,1);
% plot(n,prod(1,:),'*:');
% title('Result for x = 0');
% xlabel('n');
% ylabel('log(\Lambda (x,n)');
% 
% 
% subplot(1,3,2);
% plot(n,prod(2,:),'*:');
% title('Result for x = 0.2');
% xlabel('n');
% ylabel('log(\Lambda (x,n)');
% 
% subplot(1,3,3);
% plot(n,prod(3,:),'*:');
% title('Result for x = -1 + sqrt(2)');
% xlabel('n');
% ylabel('log(\Lambda (x,n)');



subplot(3,1,1);
plot([1 : 50 : 451],log(prod(1,:)));
legend(['t = ',num2str(0)]);
subplot(3,1,2);
plot([1 : 50 : 451],log(prod(2,:)));
legend(['t = ',num2str(0.2)]);
subplot(3,1,3);
plot([1 : 50 : 451],log(prod(3,:)));
legend(['t = ',num2str(-1 + sqrt(2))]);

figure(2)
plot([1 : 50 : 451],log(prod(1,:)),[1 : 50 : 451],log(prod(2,:)),[1 : 50 : 451],log(prod(3,:)))
hold on;
legend(['t = ',num2str(0)],['t = ',num2str(0.2)],['t = ',num2str(-1 + sqrt(2))]);

% 



%% Test for lagrangeIP and aitkeenNevilleIP
t = [0,0.5,0.7,1];
ft = [-1,1,2,1];
x = [0.25,1];
y = lagrangeIP(t,ft,x)
w = aitkenNevilleIP(t,ft,x)

%% A4 c
f = @(x) 1./(1 + x.^2);
L = [2,8,24];
I = [-3,3];
figure(1);
hold on;
for index = 1 : 3
    n = L(index);
    i = 0 : 1 : n;
    t = -3 + 6/n * i;
    ft = f(t);
    x = linspace(-3,3,200);
    lagrang = lagrangeIP(t,ft,x);
    aitken = aitkenNevilleIP(t,ft,x);
%     plot(x,lagrang);
%     plot(x,aitken);
    subplot(3,2,2 * index - 1);
    plot(x,lagrang);
    ylim([-1,1]);
    legend(['lagrange mit n = ',num2str(n)]);
    subplot(3,2,2 * index);
    plot(x,aitken);
    ylim([-1,1]);
    legend(['aitken mit n = ',num2str(n)]);
end
figure
plot(x,f(x))





%% another plot method

f = @(x) 1./(1 + x.^2);
L = [2,8,24];
I = [-3,3];
figure;
hold on;
for index = 1 : 3
    n = L(index);
    i = 0 : 1 : n;
    t = -3 + 6/n * i;
    ft = f(t);
    x = linspace(-3,3,200);
    lagrang = @(x)(lagrangeIP(t,ft,x));
    aitken = @(x) (aitkenNevilleIP(t,ft,x));
%     plot(x,lagrang);
%     plot(x,aitken);
    subplot(3,2,2 * index - 1);
    fplot(lagrang,I);
%     ylim([-1,1]);
    title(['lagrange mit n = ',num2str(n)]);
    hold on;
    fplot(f,I,'k:');
    legend('IP','f');
    subplot(3,2,2 * index);
    fplot(aitken,I);
%     ylim([-1,1]);
    hold on;
    fplot(f,I,'k:');
    title(['aitken mit n = ',num2str(n)]);
    legend('IP','f');
end



%% compute the time elapsed 
% it seems that aitken is slower.
count = 1
for  n = 16 : 50;
    i = 0 : 1 : n;
    t = -3 + 6/n * i;
    ft = f(t);
    x = linspace(-3,3,10000);
    tic
    lagrang = lagrangeIP(t,ft,x);
    temp = toc;
    time1(count) = temp;
    count = count + 1;
end

count = 1;
for  n = 16 : 50;
    i = 0 : 1 : n;
    t = -3 + 6/n * i;
    ft = f(t);
    x = linspace(-3,3,10000);
    tic;
    aitken = aitkenNevilleIP(t,ft,x);
    temp = toc;
    time2(count) = temp;
    count = count + 1;
end
figure(2)
% plot(time1,'DisplayName','t1');
plot(time1,'color','b');
hold on
% plot(time2,'DisplayName','t2');
plot(time2,'color','r');
legend('time lagrange','time aitken');