%% A5
%y(t) = exp(-1/3*t^3), y'(t) = -t^2 * y, y(0) = 1
y = @(t) exp(-1/3 * t.^3);
L = zeros(5,8);

N = 10
for n = 0 : N
    hn = 0.1 * 0.5^n;
    t = 0;
    ytj = 1;

    ytjV = 1;

    y3V = 1;
    
    y4V = 1;
   
    y5V = 1;
    while(t <= 1)
        ytj_next = ytj + hn * (-t^2) * ytj;
        ytj = ytj_next;
        
        
        ytjV_Middle = ytjV + hn/2 * (-t^2) * ytjV;
        ytjV_next = ytjV + hn * (-(t + hn/2)^2) * ytjV_Middle;
        ytjV = ytjV_next;
        
        
        % combined h/2 h/3 Euler 
        k1 = (-t^2) * y3V;
        halfV = y3V + hn/2 * k1;
        k2 = (-(t + hn/2)^2) * halfV;
        
        thirdV = y3V + hn/3 * (-t^2) * y3V;
        k3 = (-(t + hn/3)^2) * thirdV;
        
        twoThirdV = thirdV + hn/3 * (-(t + hn/3)^2) * thirdV;
        k4 = (-(t + 2/3 * hn)^2) * twoThirdV;
        y3V = y3V + hn * (-k2 + k3 + k4);
        
        % the three convergent order combined version but it seems that 
        % the order is two
        tk1 = (-t^2) * y4V;
        thalfV = y4V + hn/2 * tk1;
        tk2 = (-(t + hn/2)^2) * thalfV;
        
        tthirdV = y4V + hn/3 * (-t^2) * y4V;
        tk3 = (-(t + hn/3)^2) * tthirdV;
        
        ttwoThirdV = tthirdV + hn/3 * (-(t + hn/3)^2) * tthirdV;
        tk4 = (-(t + 2/3 * hn)^2) * ttwoThirdV;
        y4V = y4V + hn * (-2 * tk2 + 3/2 * tk3 + 3/2 * tk4);
        
        % my version . not so good :-(
        myK1 = (-t^2) * y5V;
        myK2 = (-(t+hn/2)^2)*y5V;
        myK3 = (-(t+hn/3)^2) * y5V;
        
        y5V = y5V + hn * (1/2 * myK1 + 2 * myK2 - 3/2 * myK3);
        
        
        % increase t
        t = t + hn;
        correctValue = y(t);
        diffE = abs(correctValue - ytj_next);
        diffV = abs(correctValue - ytjV_next);    
        diff23 = abs(correctValue - y3V);
        diff4 = abs(correctValue - y4V);
        diff5 = abs(correctValue - y5V);
    end
    L(1,n+1) = diffE;
    L(2,n+1) = diffV;
    L(3,n+1) =diff23;
    L(4,n+1) = diff4;
    L(5,n+1) = diff5;
end


semilogy(1:(N+1),L(1,:),'ro')
hold on
semilogy(1:(N+1),L(2,:),'bo')
semilogy(1:(N+1),L(3,:),'co')
semilogy(1:(N+1),L(4,:),'go')
semilogy(1:(N+1),L(5,:),'mo')
K1 = mean(diff(log(L(1,:)))/log(1/2));
K2 = mean(diff(log(L(2,:)))/log(1/2));
K3 = mean(diff(log(L(3,:)))/log(1/2));
K4 = mean(diff(log(L(4,:)))/log(1/2));
K5 = mean(diff(log(L(5,:)))/log(1/2));
str1 = ['euler, ord = '  num2str(K1)];
str2 = ['verbessertes euler, ord = '  num2str(K2)];
str3 = ['2,3 euler, ord = ' num2str(K3)];
str4 = ['3 order euler, ord = ' num2str(K4)];
str5 = ['my version euler, ord = ' num2str(K5)];
legend(str1,str2,str3,str4,str5)





%% A6
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
    


% %% loesung
% 
% errE = zeros(8,1)
% orderE = zeeros(8,1)
% for j = 1:8
%     
%     for i = 2: n+1
%         
%         
%     end
%     errE(j) = a