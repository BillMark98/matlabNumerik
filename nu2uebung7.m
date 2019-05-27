% % 
% % L = [0,0.2,-1 + sqrt(2)];
% % 
% % prod = zeros(3,10);
% % for index = 1 : length(L)
% %     t = L(index);
% %     count = 1;
% %     for n = 1 : 50 : 451;
% %         i = 0 : 1 : n;
% %         ti = 2 * i/n - 1;
% %         prod(index,count) = lagrangeInter(n,ti,t);
% %         count = count + 1;
% %     end
% %     
% % end  
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
% % % subplot(3,1,1);
% % plot([1 : 50 : 451],log(prod(1,:)));
% % legend(['t = ',num2str(0)]);
% % subplot(3,1,2);
% % plot([1 : 50 : 451],log(prod(2,:)));
% % legend(['t = ',num2str(0.2)]);
% % subplot(3,1,3);
% % plot([1 : 50 : 451],log(prod(3,:)));
% % legend(['t = ',num2str(-1 + sqrt(2))]);
% 
% figure(2)
% plot([1 : 50 : 451],log(prod(1,:)),[1 : 50 : 451],log(prod(2,:)),[1 : 50 : 451],log(prod(3,:)))
% hold on;
% legend(['t = ',num2str(0)],['t = ',num2str(0.2)],['t = ',num2str(-1 + sqrt(2))]);
% % calculate the i'th lagrange polynomial of degree n for the points in t
% % at point t
% function prod = lagrange_ith_Inter(i,n,ti,t)
%     prod = 1;
%     for index = 1 : (n+1)
%         if index ~= (i + 1)
%             prod = prod * (t - ti(index))/(ti(i + 1) - ti(index));
%         end
%     end
% end
% 
% % calculate the sum of the absolute value all the Lagrange poly of degree n
% % for data points ti at a specific point t
% function sum = lagrangeInter(n,ti,t)
%     sum = 0;
%     for index = 0 : n
%         sum  = sum + abs(lagrange_ith_Inter(index,n,ti,t));
%     end
% end


% 
% t = [0,0.5,0.7,1];
% ft = [-1,1,2,1];
% x = [0.25,1];
% y = lagrangeIP(t,ft,x)
% w = aitkenNevilleIP(t,ft,x)

% f = @(x) 1./(1 + x.^2);
% L = [2,8,24];
% I = [-3,3];
% figure(1);
% hold on;
% for index = 1 : 3
%     n = L(index);
%     i = 0 : 1 : n;
%     t = -3 + 6/n * i;
%     ft = f(t);
%     x = linspace(-3,3,200);
%     lagrang = lagrangeIP(t,ft,x);
%     aitken = aitkenNevilleIP(t,ft,x);
% %     plot(x,lagrang);
% %     plot(x,aitken);
%     subplot(3,2,2 * index - 1);
%     plot(x,lagrang);
%     ylim([-1,1]);
%     legend(['lagrange mit n = ',num2str(n)]);
%     subplot(3,2,2 * index);
%     plot(x,aitken);
%     ylim([-1,1]);
%     legend(['aitken mit n = ',num2str(n)]);
% end
% figure
% plot(x,f(x))





% another plot method

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



% compute the time elapsed 
% legend('lagrange n = 2','aitken n = 2','lagrange n = 8','aitken n = 8','lagrange n = 24','aitken n = 24');
% count = 1
% for  n = 16 : 22;
%     tic
%     i = 0 : 1 : n;
%     t = -3 + 6/n * i;
%     ft = f(t);
%     x = linspace(-3,3,5);
%     lagrang = lagrangeIP(t,ft,x);
% %     aitken = aitkenNevilleIP(t,ft,x);
% %     plot(x,lagrang,'DisplayName',['lagrange interpo n = ',num2str(n)]);
% %     plot(x,aitken,'DisplayName',['aitken interpo n = ',num2str(n)]);
%     subplot(3,2,2 * index - 1);
%     plot(x,lagrang,'DisplayName','lagrange interpo n = ');
% %     subplot(3,2,2 * index);
% %     plot(x,aitken,'DisplayName','aitken interpo n = ');
%     time1(count) = toc
%     count = count + 1;
% end
% count = 1;
% for  n = 16 : 22;
%     tic
%     i = 0 : 1 : n;
%     t = -3 + 6/n * i;
%     ft = f(t);
%     x = linspace(-3,3,5);
% %     lagrang = lagrangeIP(t,ft,x);
%     aitken = aitkenNevilleIP(t,ft,x);
% %     plot(x,lagrang,'DisplayName',['lagrange interpo n = ',num2str(n)]);
% %     plot(x,aitken,'DisplayName',['aitken interpo n = ',num2str(n)]);
% %     subplot(3,2,2 * index - 1);
% %     plot(x,lagrang,'DisplayName','lagrange interpo n = ');
%     subplot(3,2,2 * index);
%     plot(x,aitken,'DisplayName','aitken interpo n = ');
%     time2(count) = toc
%     count = count + 1;
% end
% figure(2)
% plot(time1,'DisplayName','t1');
% hold on
% plot(time2,'DisplayName','t2');