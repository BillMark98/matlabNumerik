
% f1 = @(y1,y2) y2;
% f2 = @(y1,y2) 8 * (1 - y1.^2) .* y2 - y1;
% f = @(y) [f1(y(1),y(2));f2(y(1),y(2))]
% N = 10;
% totalErr = zeros(1:(N + 1));
% for n = 0 : N
%     hn = 0.1 * 0.5^n;
%     y0 = 2;
%     yd0 = 0;
%     
%     y0v = [ y0;  yd0];
%     y10 = 2;
%     y20 = 0;
%     yinitv = y0v
%     count = 1;
%     
%     while(count  <=  1000)
%         y11 =  y0 + hn * 1/4 * f1(y10,y20);
%         y21 = yd0 + hn * 1/4 * f2(y10,y20);
%         
%         f_y0v = f(y0v);
%         y1v = yinitv + hn * 1/4 * f_y0v;
% %         y1v = yinitv + hn * 1/4 * f(y0v);
% %         y12 = y0 + hn * (3/32 * f1(y10,y20) + 9/32 * f1(y11,y21));
% %         y22 = yd0 + hn* (
%         f_y1v = f(y1v);
%         y2v = yinitv + hn * (3/32 * f_y0v + 9/32 * f_y1v);
% %         y2v = yinitv + hn * (3/32 * f(y0v) + 9/32 * f(y1v));
%         
%         f_y2v = f(y2v);
%         y3v = yinitv + hn * (1932/2197 * f_y0v - 7200/2197 * f_y1v + 7296/2197 * f_y2v);
% %         y3v = yinitv + hn * (1932/2197 * f(y0v) - 7200/2197 * f(y1v) + 7296/2197 * f(y2v));
% 
%         f_y3v =f(y3v);
%         y4v = yinitv + hn * (439/216 * f_y0v - 8 * f_y1v +  3680/513 * f_y2v - 845/4104 * f_y3v);
% %         y4v = yinitv + hn * (439/216 * f(y0v) - 8 * f(y1v) +  3680/513 * f(y2v) - 845/4104 * f(y3v));
%         
%         f_y4v = f(y4v);
%         ynewv = yinitv + hn * (25/216 * f_y0v + 1408/2565 * f_y2v + 2197/4104 * f_y3v - 1/5 * f_y4v);
%         
%         % for the error calculation
%         y5v = yinitv + hn * (-8/27 * f_y0v + 2 * f_y1v - 3544/2565 * f_y2v + 1859/4104 * f_y3v - 11/40 * f_y4v);
%         
%         f_y5v = f(y5v);
%         ycorv = yinitv + hn * (16/135 * f_y0v + 6656/12825 * f_y2v + 28561/56430 * f_y3v - 9/50 * f_y4v + 2/55 * f_y5v);
%         
%         error = norm(ycorv - y0v);
%         count = count + 1;
%         y0v = ynewv;
%     end
%     totalErr(n + 1) = error;
% end
%     k = mean(diff(log(totalErr))/log(1/2));
%     semilogy(1 : (N+1), totalErr,'ro');
    
%% teil 1
te = 30;
t0 = 0;
y0 = [2,0]';
err5 = zeros(5,1);
n = 350;

for i = 1 : 5
    [ya,erra] = RK45(t0,te,n,y0,@vdp8);
    err5(i) = abs(erra(n+1));
    n = 2 * n;
end
disp('the local truncation error at t = 5 for different h_i:');
err5
for i = 1 : 4
    order = log(err5(i)/err5(i+1))/log(2);
    output = ['The order of the truncation error for h',num2str(i+1), ' is ', num2str(order)]
    disp(output);
end

%% for ploting
t0 = 0;
te = 30;
n = 350;
T = linspace(t0,te, n + 1);
y0 = [2,0]';
[ya,err] = RK45(t0,te,n,y0,@vdp8);
plot(T,ya(1,:),'r.')

function dy = vdp8(t,y)
    dy = zeros(2,1);
    dy(1) = y(2);
    dy(2) = 8*(1 - y(1)^2) * y(2) - y(1);
end