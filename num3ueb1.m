% % y(t) = exp(-1/3*t^3)
% y = @(t) exp(-1/3 * t.^3);
% L = zeros(2,8);
% for n = 0 : 7
%     hn = 0.1 * 0.5^n;
%     t = 0;
%     ytj = 1;
%     err = 0;
%     
%     ytjV = 1;
%     errV = 0;
%     while(t <= 0.1)
%         ytj_next = ytj + hn * (-t^2) * ytj;
%         ytj = ytj_next;
%         
%         
%         ytjV_Middle = ytjV + hn/2 * (-t^2) * ytjV;
%         ytjV_next = ytjV + hn * (-(t + hn/2)^2) * ytjV_Middle;
%         ytjV = ytjV_next;
%         t = t + hn;
%         correctValue = y(t);
%         diffE = abs(correctValue - ytj_next);
%         if(diffE > err) 
%             err = diffE;
%         end
%         diffV = abs(correctValue - ytjV_next);
%         if(diffV > errV)
%             errV = diffV;
%         end
%     end
%     L(1,n+1) = err;
%     L(2,n+1) = errV;
% end
% 
% semilogy(1:8,L(1,:),'ro')
% hold on
% semilogy(1:8,L(2,:),'bo')
% 
% k1 = mean(diff(log(L(1,:)))/log(1/2));
% k2 = mean(diff(log(L(2,:)))/log(1/2));
% str1 = ['euler, ord = '  num2str(k1)];
% str2 = ['verbessertes euler, ord = '  num2str(k2)];
% legend(str1,str2)





% A6
% [x,y]' = 
dist = 1;
v0 = [1,0]';
h = 0.1 * 0.5^7;
t = 0;
f = @(t,x,y) 2/norm([t - y,x]) * [-x,t-y]';
v = zeros(2,1);
while(dist > 1e-5)
    % choose the x-coord as relative dist
    % since the norm of the distance will be inaccurate
    % by examining the vector returned  [0,0.6646] we see that
    % x-coord already 0, while y-coord and t has a large difference
    % but norm([0,t] - v) in this case is exactly abs(t - v(2)) and 
    % cannot be smaller than 1e-5 
    % so the function will run forever. Therefore, use x-coord as dist
    % instead
    v = v0 + h/2 * f(t,v0(1),v0(2));
    v = v0 + h * f(t+h/2, v(1),v(2));
    t = t + h;
    v0 = v;
    dist = v(1);
end
v
    