%% first equation
format long
A = [-2,1;1,-2];
g = @(t) [2 * sin(t);2 *(cos(t) - sin(t))];
y0 = [2,3];
t0 = 0;
tend = 10;
y_ans = @(t) 2 * exp(-t) .* [1;1] + [sin(t);cos(t)];
correctAns = y_ans(10);
% mode = 1
mode = 1;
epsilon = 1e-1;
flagFound = false;
for n = 50:100
%     try 
%         y = linearSDIRK(mode, A,g,n,y0,t0,tend);
%     catch ME
%         continue;
%     end
    y = linearSDIRK(mode, A,g,n,y0,t0,tend);
    if(abs(y(:,n + 1) - correctAns) <= epsilon)
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
    h1 = h + 0.1;
    n1 = round((tend - t0)/h1);
    T2 = linspace(t0,tend,n1 + 1);
    yh1 = linearSDIRK(mode, A,g,n1,y0,t0,tend);
    subplot(2,2,1)
    plot(T,y(1,:),'b-');
    title('the first component of the answer with h');
    subplot(2,2,2)
    plot(T2,yh1(1,:),'r');
    title('the first component of the answer with larger h');
    subplot(2,2,3)
    plot(T,y(2,:),'b-');
    title('the second component of the answer');
    subplot(2,2,4)
    plot(T2,yh1(2,:),'r');
    title('the second component of the answer with larger h');
end
% mode = 2
    

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
A = [-2,1;998,-999];
g = @(t) [2 * sin(t);999 * (cos(t) - sin(t))];
y0 = [2,3];
t0 = 0;
tend = 10;
y_ans = @(t) 2 * exp(-t) .* [1;1] + [sin(t);cos(t)];
correctAns = y_ans(10);
% mode = 1
mode = 1;
epsilon = 1e-3;
flagFound = false;
for n = 100:140
%     try 
%         y = linearSDIRK(mode, A,g,n,y0,t0,tend);
%     catch ME
%         continue;
%     end
    y = linearSDIRK(mode, A,g,n,y0,t0,tend);
    if(abs(y(:,n + 1) - correctAns) <= epsilon)
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
    h1 = h + 0.1;
    n1 = round((tend - t0)/h1);
    T2 = linspace(t0,tend,n1 + 1);
    yh1 = linearSDIRK(mode, A,g,n1,y0,t0,tend);
    subplot(2,2,1)
    plot(T,y(1,:),'b-');
    title('the first component of the answer with h');
    subplot(2,2,2)
    plot(T2,yh1(1,:),'r');
    title('the first component of the answer with larger h');
    subplot(2,2,3)
    plot(T,y(2,:),'b-');
    title('the second component of the answer');
    subplot(2,2,4)
    plot(T2,yh1(2,:),'r');
    title('the second component of the answer with larger h');
end