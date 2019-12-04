%% first equation
format long
A = [-2,1;1,-2];
g = @(t) [2 * sin(t);2 *(cos(t) - sin(t))];
y0 = [2,3];
t0 = 0;
tend = 10;
y_ans = @(t) 2 * exp(-t) .* [1;1] + [sin(t);cos(t)];
epsilon = 1e-3;
% correctAns = y_ans(10);
% here view the calculated answer as the exact solution if 
% the absolute error at each time is less than epsilon
% i.e max(abs(y - y_ans)) <= epsilon ( use maximum norm)
% since in this case, the eigenvalue of A is rather small
% so with small variation of stepsize, the h * lambda as the notation in the 
% lecture, still lie within the stable interval, so the answer with larger
% step size is still close to the exact solution
for mode = 1 : 2
    if(mode == 1)
        BETA = "1/2 + $$\sqrt{3}$$/6";
    else
        BETA = "1/2 - $$\sqrt{3}$$/6";
    end
    
    flagFound = false;
    for n = 50:1000
    %     try 
    %         y = linearSDIRK(mode, A,g,n,y0,t0,tend);
    %     catch ME
    %         continue;
    %     end
        y = linearSDIRK(mode, A,g,n,y0,t0,tend);
        T = linspace(t0,tend,n + 1);
        if(max(abs(y - y_ans(T))) <= epsilon)
            disp('minimum h found: ');
            h = (tend - t0)/n
            disp(h);
            flagFound = true;
            break;
        end
    end
    if(~flagFound)
        error('step need to be smaller');
    else
        figure;
        T = linspace(t0,tend,n + 1);
        h1 = h + 0.01;
%         h1 = (tend - t0)/n1;
        n1 = round((tend - t0)/h1);
        T2 = linspace(t0,tend,n1 + 1);
        yh1 = linearSDIRK(mode, A,g,n1,y0,t0,tend);
        subplot(2,1,1)
        plot(T,y(1,:),'b-');
        title('the first component of the answer');
        hold on
        plot(T2,yh1(1,:),'r');
        legend(['correct y_1, step = ',num2str(h)],['y_1 with larger step = ',num2str(h1)])
        subplot(2,1,2)
        plot(T,y(2,:),'b-');
        title('the second component of the answer');
        hold on
        plot(T2,yh1(2,:),'r');
        legend(['correct y_2, step = ',num2str(h)],['y_2 with larger step = ',num2str(h1)],'Position',[0.3,0.3,0.2,0.1])
        myTitle = sgtitle(["SDIRK-method with beta = ", BETA]);
        set(myTitle,'Interpreter','latex','fontsize',12);
    end
end

% 
%     figure;
%     T = linspace(t0,tend,n + 1);
%     yAns = y_ans(T);
%     y1 = yAns(1,:);
%     y2 = yAns(2,:);
%     subplot(2,1,1)
%     plot(T,y1,'b-');
%     hold on
%     plot(T,y(1,:),'r');
%     legend('correct y_1','calculated y_1')
%     title('the first component of the answer');
%     subplot(2,1,2)
%     plot(T,y2,'b-');
%     hold on
%     plot(T,y(2,:),'r');
%     legend('correct y_2','calculated y_2');
%     title('the second component of the answer');





% figure;
%     T = linspace(t0,tend,n + 1);
%     yAns = y_ans(T);
%     y1 = yAns(1,:);
%     y2 = yAns(2,:);
%     subplot(2,2,1)
%     plot(T,y1,'b-');
%     title('the correct first component of the answer');
%     subplot(2,2,2)
%     plot(T,y(1,:),'r');
%     title('the calculated first component of the answer');
%     subplot(2,2,3)
%     plot(T,y2,'b-');
%     title('the correct second component of the answer');
%     subplot(2,2,4)
%     plot(T,y(2,:),'r');
%     title('the calculated second component of the answer');

%% second equation
% we see that when we change step a little bit when beta = 1/2 - sqrt(3)/6,
% the answer varies a lot. It's because we have calculated in the a) part
% that the stable interval for beta = 1/2 - sqrt(3)/6 is smaller than the
% intervall of that of beta = 1/2 + sqrt(3)/6.
% in partitcular the larger step size h1 * max(abs(eig(A))) lies outside
% the stable interval, so the answer will explode. ( see the calculation at
% the end). Note that the stable interval for beta = 1/2 - sqrt(3)/6 is
% (-(4 * sqrt(3) + 6),0)
A = [-2,1;998,-999];
g = @(t) [2 * sin(t);999 * (cos(t) - sin(t))];
y0 = [2,3];
t0 = 0;
tend = 10;
y_ans = @(t) 2 * exp(-t) .* [1;1] + [sin(t);cos(t)];
correctAns = y_ans(10);
epsilon = 1e-3;
flagFound = false;
for mode = 1 : 2
    if(mode == 1)
        BETA = "1/2 + $$\sqrt{3}$$/6";
    else
        BETA = "1/2 - $$\sqrt{3}$$/6";
    end
    epsilon = 1e-3;
    flagFound = false;
    for n = 50:1000
    %     try 
    %         y = linearSDIRK(mode, A,g,n,y0,t0,tend);
    %     catch ME
    %         continue;
    %     end
        y = linearSDIRK(mode, A,g,n,y0,t0,tend);
        T = linspace(t0,tend,n + 1);
        if(max(abs(y - y_ans(T))) <= epsilon)
            disp('minimum h found: ');
            h = (tend - t0)/n
            disp(h);
            flagFound = true;
            break;
        end
    end
    if(~flagFound)
        error('step need to be smaller');
    else
        figure;
        T = linspace(t0,tend,n + 1);
%         n1 = n - 30;
%         h1 = (tend - t0)/n1;
        h1 = min(h + 0.0005, (4 * sqrt(3) + 7)/1000);
%         h1 = (tend - t0)/n1;
        n1 = round((tend - t0)/h1);
        T2 = linspace(t0,tend,n1 + 1);
        yh1 = linearSDIRK(mode, A,g,n1,y0,t0,tend);
        subplot(2,1,1)
        plot(T,y(1,:),'b-');
        title('the first component of the answer');
        hold on
        plot(T2,yh1(1,:),'r');
        legend(['correct y_1, step = ',num2str(h)],['y_1 with larger step = ',num2str(h1)])
        subplot(2,1,2)
        plot(T,y(2,:),'b-');
        title('the second component of the answer');
        hold on
        plot(T2,yh1(2,:),'r');
        legend(['correct y_2, step = ',num2str(h)],['y_2 with larger step = ',num2str(h1)],'Position',[0.3,0.3,0.2,0.1])
        myTitle = sgtitle(["SDIRK-method with beta = ", BETA]);
        set(myTitle,'Interpreter','latex','fontsize',12);
    end
end


% for n = 500:1000
% %     try 
% %         y = linearSDIRK(mode, A,g,n,y0,t0,tend);
% %     catch ME
% %         continue;
% %     end
%     y = linearSDIRK(mode, A,g,n,y0,t0,tend);
%     if(abs(y(:,n + 1) - correctAns) <= epsilon)
%         disp('minimum h found: ');
%         h = (tend - t0)/n
%         disp(h);
%         flagFound = true;
%         break;
%     end
% end
% if(~flagFound)
%     error('step need to be smaller');
% else
%     figure;
%     T = linspace(t0,tend,n + 1);
%     h1 = h + 0.00044;
%     n1 = round((tend - t0)/h1);
%     T2 = linspace(t0,tend,n1 + 1);
%     yh1 = linearSDIRK(mode, A,g,n1,y0,t0,tend);
%     subplot(2,2,1)
%     plot(T,y(1,:),'b-');
%     title('the first component of the answer with h');
%     subplot(2,2,2)
%     plot(T2,yh1(1,:),'r');
%     title('the first component of the answer with larger h');
%     subplot(2,2,3)
%     plot(T,y(2,:),'b-');
%     title('the second component of the answer');
%     subplot(2,2,4)
%     plot(T2,yh1(2,:),'r');
%     title('the second component of the answer with larger h');
% end


h1 * max(abs(eig(A))) - (4 * sqrt(3) + 6)
h * max(abs(eig(A))) - (4 * sqrt(3) + 6)