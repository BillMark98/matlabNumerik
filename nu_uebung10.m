% example from the help document
% remember that the coeff are local 
% which means given interval [x1,x2], and coeff [a,b,c,d]
% the polynomial is given by a(x-x1)^3 + b(x-x1)^2 + c(x-x1) + d
% so in this case [-8,-4] cc = [-1/4,1,0]
% the polynomial is -1/4 (x + 8)^2 + (x + 8)

% subplot(2,2,1)
% cc = [-1/4 1 0]; 
% pp1 = mkpp([-8 -4],cc);
% xx1 = -8:0.1:-4; 
% plot(xx1,ppval(pp1,xx1),'k-')
% 
% subplot(2,2,2)
% pp2 = mkpp([-4 0],-cc);
% xx2 = -4:0.1:0; 
% plot(xx2,ppval(pp2,xx2),'k-')
% 
% subplot(2,1,2)
% pp = mkpp([-8 -4 0 4 8],[cc;-cc;cc;-cc]);
% xx = -8:0.1:8;
% plot(xx,ppval(pp,xx),'k-')
% hold on
% line([-4 -4],ylim,'LineStyle','--')
% line([0 0],ylim,'LineStyle','--')
% line([4 4],ylim,'LineStyle','--')
% hold off


cc = [1 0 1 0 -1];
x0 = -1;
x1 = 0;
x2 = 1;
pp = mkpp([x0,x1,x2],[cc;-cc]);
xx = x0:0.1:x2;
plot(xx,ppval(pp,xx),'k-')
hold on
line([0,0],ylim,'LineStyle','--')
hold off


% f = [1;1;2;2];
% a = 0;
% b = 3;
% pp = cubicSpline(f,a,b);
% xx = a:0.1:b;
% subplot(2,1,1);
% plot(xx,ppval(pp,xx),'k-')
% hold on
% ylim([0,2.5]);
% line([1,1],ylim,'LineStyle','--')
% line([2,2],ylim,'LineStyle','--')
% line([3,3],ylim,'LineStyle','--')
% 
% subplot(2,1,2)
% x = 0 : 3
% pps = spline(x,f)
% plot(xx,ppval(pps,xx),'k-')
% ylim([0,2.5]);
% hold on
% line([1,1],ylim,'LineStyle','--')
% line([2,2],ylim,'LineStyle','--')
% line([3,3],ylim,'LineStyle','--')

% f = @(x) 1./(1 + x.^2);
% n = 11;
% x = linspace(-5,5,n);
% I = [-5,5]
% ff = f(x);
% a = -5;
% b = 5;
% pp = cubicSpline(ff,a,b);
% xx = a : 0.1 : b;
% plot(xx,ppval(pp,xx),'r-')
% hold on
% fplot(f,I,'k:')
% legend('spline','f')