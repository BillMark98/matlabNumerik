% 
% L = [0,0.2,-1 + sqrt(2)];
% 
% prod = zeros(3,10);
% for index = 1 : length(L)
%     t = L(index);
%     count = 1;
%     for n = 1 : 50 : 451;
%         i = 0 : 1 : n;
%         ti = 2 * i/n - 1;
%         prod(index,count) = lagrangeInter(n,ti,t);
%         count = count + 1;
%     end
%     
% end  
% 
% 
% % subplot(3,1,1);
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



t = [0,0.5,0.7,1];
ft = [-1,1,2,1];
x = [0.25,1];
y = lagrangeIP(t,ft,x)
w = aitkenNevilleIP(t,ft,x)

f = @(x) 1./(1 + x.^2);
L = [2,8,24];
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
    plot(x,lagrang);
    plot(x,aitken);
%     subplot(3,2,2 * index - 1);
%     plot(x,lagrang);
%     subplot(3,2,2 * index);
%     plot(x,aitken);
end
legend('lagrange n = 2','aitken n = 2','lagrange n = 8','aitken n = 8','lagrange n = 24','aitken n = 24');
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