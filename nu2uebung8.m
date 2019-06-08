f = @(x) 1./(1 + x.^2);
L = [5,10,40];
ende = 2;
I = [-ende,ende];
figure;

for index = 1 : 3
    n = L(index);
    left = -ende;
    subplot(1,3,index);
    
    fplot(f,I,'k:');
    hold on;
    i = 0 : 1 : n;
    x = -ende + 2 * ende *i/n;
    y = f(x);
    p = polyfit(x,y,n);
    interpol = @(x) polyval(p,x);
    fplot(interpol,I);
%     ylim([-5,1.5]);
    title(['equidist IP n = ',num2str(n)],'Units','normalized','Position',[0.68,1,0]);
    legend('f','IP');
end

format long;
figure;
for index = 1 : 3
    n = L(index);
    subplot(1,3,index);
    fplot(f,I,'k:');
    hold on;
    i = 0 : 1 : n;
    x = ende * cos((2*i+1)/(2*n+2)*pi);
    y = f(x);
    p = polyfit(x,y,n);
    interpol = @(x) polyval(p,x);
    fplot(interpol,I);
%     ylim([-5,1.5]);
    title(['tscheby points IP n = ',num2str(n)]);
    legend('f','IP');
end

figure;
for index = 1 : 3
    n = L(index);
    subplot(1,3,index);
    hold on;
    i = 0 : 1 : n;
    
    x1 = -ende + 2 * ende *i/n;
    p1 = poly(x1);
    x2 = ende * cos((2*i+1)/(2*n+2)*pi);
    p2 = poly(x2);
  
    knotenpoly1 = @(x) polyval(p1,x);
    fplot(knotenpoly1,I);
%     ylim([-5,1.5]);
    knotenpoly2 = @(x) polyval(p2,x);
    fplot(knotenpoly2,I);
    t = title(['knotenpoly n = ',num2str(n)],'Units','normalized','Position',[0.68,1.01,0]);
    
    legend('equidist','tscheby');
end