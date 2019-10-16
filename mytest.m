% % delta0 = zeros(9,1);
% % delta1 = zeros(9,1);
% % count = 0;
% % for n = 10:10:100
% %     count = count + 1;
% %     A = eye(n);
% %     A(:,n) = 1;
% %     for i = 1 : n
%         for j = 1 : (i-1)
% %             A(i,j) = -1;
% %         end
% %     end
% %     b = zeros(n,1);
% %     result = zeros(n,1);
% %     for i = 1 : n
% %         b(i) = (i-1)/n;
% %         result(i) = -1/n;
% %     end
% %     b(n) = 1;
% %     result(n) = 1/n;
% % 
% %     [LR,p] = LR_Decomp(A);
% %     x = LR_Solve(LR,p,b);
% %     x1 = lr_nachiteration(A,LR,p,b,10);
% %     delta0(count) = norm(x-result)/norm(result);
% %     delta1(count) = norm(x1-result)/norm(result);
% % end
% % N = 10:10:100;
% % figure
% % semilogy(N,delta0);
% % hold on
% % semilogy(N,delta1);
% % legend('delta0','delta1');
% % 
% % n = 3;
% % p = 1:n;
% % M = n*mod(p'+p-(n+3)/2,n) + mod(p'+2*p-2,n) + 1;
% 
% 
% 
% % bg = [0 0 85] ;
% % fg = [255 255 255];
% % sz = get(0,'screensize');
% % rand('state',0)
% % X = finitefern(500000,sz(4),sz(3));
% % d = fg - bg;
% % R = uint8(bg(1) + d(1)*X);
% % G = uint(bg(2) + d(2)*X);
% % B = uint(bg(3) + d(3)*X);
% % F = cat(3,R,G,B);
% % % imwrite(F,'myfern.png','png','bitdepth',8);
% % 
% % 
% % A = [21,-51,86;0,-25,62;28,32,-77]
% % res = houseQR(A)
% % G1 = [3/5 0 4/5; 0 1 0;-4/5 0 3/5]
% % G2 = [1 0 0; 0 -5/13 12/13; 0 -12/13 -5/13]
% % G2*G1*A
% % 
% % format long;
% % [K,T] = arrhenius_data();
% % R = 8.314;
% % E1 = 2.75e5;
% % A1 = arrhenius_A(K,T,E1);
% % [A2,E2] = arrhenius_AE(K,T);
% % 
% % plot(T,K,'.');
% % hold on;
% % f1 = arrayfun(@(x) A1*exp(-E1/(R*x)),T);
% % f2 = arrayfun(@(x) A2*exp(-E2/(R*x)),T);
% % plot(T,f1);
% % plot(T,f2);
% % xlabel('T');
% % ylabel('K');
% % legend('data','f1','f2');
% % figure(2);
% % semilogy(T,K,'.');
% % hold on;
% % semilogy(T,f1);
% % semilogy(T,f2);
% % xlabel('T');
% % ylabel('K');
% % legend('data','f1','f2');
% % figure(3);
% % T1 = arrayfun(@(x) 1000/x,T);
% % semilogy(T1,K,'.');
% % hold on;
% % semilogy(T1,f1);
% % semilogy(T1,f2);
% % xlabel('1000/T');
% % ylabel('K');
% % legend('data','f1','f2');
% % 
% % norm(K-f1)^2
% % norm(K-f2)^2
% 
% 
% % x = 0 : 1 : 10;
% % y = 0 : 0.1 : 10;
% % z1 = sin(2*3*pi/5*x);
% % z2 = sin(2*3*pi/5*y);
% % plot(x,z1,'o');
% % hold on
% % plot(y,z2);
% % %
% % z3 = sin(-2*2*pi/5*x);
% % hold on;
% % plot(x,z3,'*');
% % plot(y,z4);
% 
% % n = 4
% % A = magic(n)
% % p = randperm(n)
% % q = randperm(n)
% % A = A(p,q)
% % sum(A)
% % sum(A')'
% % sum(diag(A))
% % sum(diag(flipud(A)))
% % rank(A)
% 
% 
% % A = [3,7;0,12;4,1]
% % b = [10;1;5]
% % QR = houseQR(A)
% % [x,res] = houseQR_Solve(QR,b)
% 
% % 
% % delta = 1e-6;
% % A = [sqrt(2),sqrt(2);sqrt(2)+delta, sqrt(2);delta,0;0,delta]
% % b = [2*sqrt(2);2*sqrt(2)+delta;delta;delta]
% % QR = houseQR(A)
% % [x,res] = houseQR_Solve(QR,b)
% % norm(A*x-b)
% % x1 = A\b
% % 
% % Loes = [1;1]
% % norm(Loes - x)/norm(Loes)
% % x2 = (A'*A)\(A'*b)
% % norm(Loes - x2)/norm(Loes)
% 
% 
% % N = 7;
% % N1 = 2;
% % index = floor(N/2);
% % k = -index : index;
% % A = k;
% % for i = 1 : N
% %     if k(i) ~= 0
% %         A(i) = 1/N*(sin(2*pi*k(i)*(N1 + 1/2)/N)/sin(pi*k(i)/N))^2;
% %     end
% % end
% % A(index + 1) = (2*N1+1)^2/N;
% % x = -index : index;
% % y = x;
% % for j = 1 : N
% %     sum = 0;
% %     for i = (index + 2) : N
% %         sum = sum + A(i) * 2 * cos(2*pi/N*k(i)*x(j));
% %     end
% %     sum = sum + A(index + 1);
% %     y(j) = sum;
% % end
% % 
% % k = 1;
% % X = (-index - k * N) : (index + k * N);
% % Y = repmat(y,1,(2*k+1));
% % 
% % 
% % stem(X,Y);
% 
% % a = 0.7
% % bb = [1, -a]; %-- Filter Coefficients
% % ww = -2*pi:(pi/100):2*pi; %-- omega hat
% % HH = freqz(1,bb, ww); %<--freekz.m is an alternative
% % subplot(2,1,1);
% % plot(ww, abs(HH))
% % subplot(2,1,2);
% % plot(ww, angle(HH))
% % xlabel('Normalized Radian Frequency')
% 
% % A = [2 1;1 2]
% % theta = [0:2*pi/50:2*pi];
% % circle = [cos(theta);sin(theta)];
% % ellipse = A*circle;
% % axis([-4 4 -4 4]); axis('square')
% % plot(circle(1,:),circle(2,:),ellipse(1,:),ellipse(2,:))
% 
% % figure(2)
% % H = [ -6 -6 -7 0 7 6 6 -3 -3 0 0 -6;-7 2 1 8 1 2 -7 -7 -2 -2 -7 -7]
% % A = [0,1;1,0];
% % HH = A * H;
% % x = HH(1,:)'; y = HH(2,:)';
% % axis([-10 10 -10 10]);axis('square')
% % plot(x,y,'o',x,y,'-')
% 
% 
% 
% % alpha = -1/3;miu = 2 * 4*pi*1e-7;epsilon = 8*8.854*1e-12;
% % omega = 1e08;
% % temp = ((alpha/omega)^2/(miu*epsilon/2)+1)
% % beta = omega*sqrt((miu*epsilon/2*(temp + 1)))
% % thet = sqrt(temp^2 - 1)
% % sqrt(miu/epsilon)/sqrt(temp)
% % atan(thet)/2 * 180 / pi
% 
% 
% % N = 10;
% % Res = zeros(N + 1,5);
% % a0 = 2.5;
% % a_prev = a0;
% % a_pp = a_prev;
% % a = sqrt(5);
% % Res(1,1) = a0;
% % Res(1,5) = a - a0;
% % for i = 2 :  (N + 1)
% % %     a_next = 1 + a_prev - 1/5 * a_prev^2;
% %     a_next = 1/2*(a_prev + 5/a_prev);
% %     Res(i,1) = a_next;
% %     if i > 2 
% %         Res(i,2) = (a_next - a_prev)/(a_prev - Res(i-2,1));
% %         Res(i,3) = Res(i,2)/(1 - Res(i,2)) * (a_next - a_prev);
% %     end
% %     Res(i - 1, 4) = a_next - a_prev;
% %     Res(i,5) = a - a_next;
% %     a_prev = a_next;
% % end
% % Res
% 
% % a0 = 1.2;
% % a1 = atan(2*a0);
% % a_prev= a0;
% % a = a1;
% % L = 2/5;
% % err = L/(1 - L) * (a - a_prev);
% % while err > 1e-4
% %     a_prev = a;
% %     a = atan(2 * a_prev);
% %     err = L/(1 - L) * (a - a_prev);
% %     disp(a);
% %     disp(err);
% % end
% % 
% % x0 = 4;
% % x_prev = x0;
% % f = @(x) cos(x) + exp((x - pi)^2);
% % df = @(x) -sin(x) + 2*(x - pi)*exp((x - pi)^2);
% % for i = 1 : 3
% %     
% %     x_ne = x_prev - 2 * f(x_prev)/df(x_prev);
% %     x_prev = x_ne
% %     f(x_ne)
% % end
% 
% % delta_t = 0.15;
% % f1 = @(t) t;
% % f2 = @(t) exp(t);
% % f3 = @(t) t^3;
% % f4 = @(t) sin(t);
% % A = zeros(6,3);
% % AHat = zeros(6,4);
% % x = [1.2,0.6,0.6]'
% % xHat = [1.2,0.6,1.6,0.9]'
% % b = zeros(6,1);
% % bHat = zeros(6,1);
% % for i = 1 : 6
% %     A(i,:) = [f1(delta_t*i),f2(delta_t*i),f3(delta_t*i)];
% %     AHat(i,:) = [f1(delta_t*i),f2(delta_t*i),f3(delta_t*i),f4(delta_t*i)];
% %     b(i,1) = x(1)*f1(delta_t * i) + x(2) * f2(delta_t * i) + x(3) * f3(delta_t * i);
% %     bHat(i,1) = xHat(1)*f1(delta_t * i) + xHat(2) * f2(delta_t * i) + xHat(3) * f3(delta_t * i) + xHat(4) * f4(delta_t * i);
% % end
% % A
% % AHat
% % A*x - b
% % AHat * xHat - bHat
% % cond(A)
% % cond(AHat)
% % 
% % delta_b1 = 1e-2 * [1,0,-1,-1,-0.5,1]'
% % delta_b2 = 1e-2 * [-1,1,1,-0.5,-2,1]'
% % b = b + delta_b1
% % bHat = bHat + delta_b2
% % xnew = A\b
% % xnew = (A'*A)\(A'*b)
% % 
% % xHatnew = (AHat'* AHat)\(AHat' * bHat)
% % rel1 = norm(x - xnew)/norm(x)
% % rel2 = norm(xHat - xHatnew)/norm(xHat)
% % 
% % [U,S,V] = svd(A)
% % [UHat,SHat,VHat] = svd(AHat)
% % 
% % b_cord = U\delta_b1
% % bHat_cord = UHat\delta_b2
% 
% 
% 
% x = 0: 0.01 : 3;
% % y = x.^4 - x -1;
% % z = x * 0;
% % plot(x,y,x,z)
% f1 = @(x) x^4 - 1;
% f2 = @(x) (x + 1)^(1/4);
% f3 = @(x) (3*x^4 + 1)/(4*x^3 - 1);
% % x0 = 0.9;
% % y1 = zeros(1,10);
% % y2 = zeros(1,10);
% % y3 = zeros(1,10);
% % y1(1) = x0;
% % y2(1) = x0;
% % y3(1) = x0;
% % for i = 1 : 9
% %     y1(i + 1) = f1(y1(i)); 
% %     y2(i + 1) = f2(y2(i)); 
% %     y3(i + 1) = f3(y3(i)); 
% % end
% % y1
% % y2
% % y3
% 
% % syms F2(x) F3(x)
% % F2(x) = (x + 1)^(1/4);
% % df2 = diff(F2,x);
% % F3(x) = (3*x^4 + 1)/(4*x^3 - 1);
% % df3 = diff(F3,x);
% % L2 = abs(eval(df2(1.2207)))
% % L22 = abs(eval(df2(0.9)))
% % L3 = abs(eval(df3(1.2207)))
% % L33 = abs(eval(df3(0.9)))
% % y2_Ak = zeros(1,10);
% % y2_diff = zeros(1,10);
% % y3_Ak = zeros(1,10);
% % y3_diff = zeros(1,10);
% % 
% % for i = 3 : 10
% %     
% %     y2_Ak(i) = (y2(i)-y2(i-1)) / (y2(i-1) - y2(i-2));
% %     y2_diff(i) = y2_Ak(i)/(1 - y2_Ak(i)) * (y2(i) - y2(i-1));
% %     y3_Ak(i) = (y3(i) - y3(i-1)) / (y3(i-1) - y3(i-2));
% %     y3_diff(i) = y3_Ak(i)/(1 - y3_Ak(i)) * (y3(i) - y3(i-1));
% %     
% % end
% % semilogy([3:10],y2_Ak(3:10))
% % semilogy([3:10],y3_Ak(3:10))
% 
% % f = @(x) cos(x)*cosh(x) + 1
% % x0 = 2.25;
% % f1 = @(x) -sin(x)*cosh(x) + cos(x)*sinh(x)
% % x1 = zeros(1,5);
% % x2 = zeros(1,5);
% % x1(1) = x0;
% % x2(1) = x0;
% % for i = 2 : 5
% %     x1(i) = x1(i-1) - f(x1(i-1))/f1(x1(i-1));
% %     x2(i) = x2(i-1) - f(x2(i-1))/f1(x2(i-1)) - 1/(2*(f1(x2(i-1)))^2) * (f(x2(i-1)))^2;
% % end
% % x1
% % x2
% 
% % x0 = [1;1]
% % F = @(x1,x2) [f1(x1,x2);f2(x1,x2)]
% % b = F(x0(1),x0(2))
% % dF = diff(F,[x1,x2])
% 
% 
% 
% % 
% % %CN01.m
% % 
% % % V/I connections for a resistor, capacitor and inductor
% % % Finite Difference Method for circuit anlysis
% % 
% % % Ian Cooper
% % % email: ian.cooper@sydney.edu.au
% % % School of Physics, University of Sydney
% % 
% % % DOING PHYSICS WITH MATLAB 
% % %    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% % % 180121
% % 
% % clear all
% % close all
% % clc
% % 
% % % INPUTS: S.I. UNITS ====================================================
% %   R = 1.0e3;
% %   C = 1e-6;
% %   L = 10e-3;
% %   IR = 0.10;
% %   IC = 0.10;
% %   VL = 10;
% %   f = 1e3;
% %   N = 200;
% % 
% % % CALCULATIONS ==========================================================  
% %   T = 1/f;
% %   w = 2*pi*f;
% %   tMax = 2*T;
% %   t = linspace(0,tMax,N);
% %   dt = t(2)-t(1);
% %   
% %   iR = IR.*cos(w*t);
% %   vR = R .* iR;
% %   
% %   kC = 2*dt/C;
% %   iC = IC.* cos(w*t);
% %   vC = zeros(1,N);
% %   vC(2) = vC(1)+(dt/C)*iC(1);
% %   
% %   kL = 2*dt/L;
% %   vL = VL .* sin(w*t);
% %   iL = zeros(1,N);
% %   iL(2) = iL(1)+(dt/L)*vL(1);
% % 
% % for c = 3 : N-1
% %  vC(c) = vC(c-2) + kC*iC(c-1);
% %  iL(c) = iL(c-2) + kL*vL(c-1);
% %  %iL(c) = iL(c-1)+(dt/L)*vL(c-1);
% % end
% %  
% %  iL = iL - max(iL)/2;   % not sure why need to shift the current
% %  
% % % Peak values /Impedances / reactances
% %   VR = max(vR);
% %   IR = max(iR);
% %   VC = max(vC);
% %   IC = max(iC);
% %   VL = max(vL);
% %   IL = max(iL);
% % 
% %   ZR = VR/IR;        % resistance
% %   ZC = VC/IC;        % capacitive reactance
% %   ZL = VL/IL;        % inductive reactance
% % 
% % % Theorteical values for impedance
% %   ZR_T = VR/IR;
% %   ZC_T = 1/(w*C);
% %   ZL_T = w*L;
% % 
% % 
% % % =====================================================================
% % %   DISPLAY RESULTS
% % % =====================================================================
% %      fprintf(' input frequency    f =  %3.2f   Hz \n',f);
% %      disp('  ') 
% % disp(' RESISTOR  ');
% %      fprintf(' resistance   R =  %3.2f   ohms \n',R);
% %      fprintf(' peak voltage VR =  %3.2f  V  \n',VR);
% %      fprintf(' peak current IR =  %3.2f  mA  \n',1e3*IR);
% %      fprintf(' impedance    ZR =  %3.2f  ohms  \n',ZR);
% %      fprintf(' theoretical impedance ZR_T =  %3.2f ohms \n',ZR_T);
% %      disp('  ') 
% % disp('  CAPACITOR  ');
% %      fprintf(' capacitance  C =  %3.2e  F \n',C);
% %      fprintf(' peak voltage VC =  %3.2f  V  \n',VC);
% %      fprintf(' peak current IC =  %3.2f  mA  \n',1e3*IC);
% %      fprintf(' impedance    ZC =  %3.2f  ohms  \n',ZC);
% %      fprintf(' theoretical impedance ZC_T =  %3.2f ohms \n',ZC_T);
% %      disp('  ') 
% % disp('  INDUCTOR  ');
% %      fprintf(' inductance   L =  %3.2e  H \n',L);
% %      fprintf(' peak voltage VL =  %3.2f  V  \n',VL);
% %      fprintf(' peak current IL =  %3.2f  mA  \n',1e3*IL);
% %      fprintf(' impedance    ZL =  %3.2f  ohms  \n',ZL);
% %      fprintf(' theoretical impedance ZL_T =  %3.2f ohms \n',ZL_T);
% %      
% % % ======================================================================
% % % GRAPHICS
% % % ======================================================================
% % 
% % figure(1)    % RESISTOR ------------------------------------------------
% %    FS = 14;
% %    pos = [0.07 0.05 0.28 0.32];
% %    set(gcf,'Units','normalized');
% %    set(gcf,'Position',pos);
% %    set(gcf,'color','w'); 
% %       
% %    yyaxis left 
% %    xP = 1e3.*t; yP = vR;
% %    plot(xP,yP,'b','linewidth',2');
% %    grid on
% %    xlabel('time  t  [ms]');
% %    ylabel('v_R  [ V ]');
% %    
% %    
% %    ax = gca;
% %    ax.YAxis(1).Color = 'b';
% %    set(gca,'fontsize',FS);
% %    set(gca,'xLim',[0 1e3*max(t)]);
% %    
% %    yyaxis right
% %    xP = 1e3.*t; yP = iR.*1e3;
% %    plot(xP,yP,'r--','linewidth',2');
% %    ylabel('i_R  [ mA ]');
% %    ax = gca;
% %    ax.YAxis(2).Color = 'r';
% %    
% %    tm  = ' RESISTOR: voltage and current in phase'; 
% %    title(tm,'fontweight','normal','fontsize',12);
% % 
% % 
% % figure(2);  % CAPACITOR -------------------------------------------------
% %    FS = 14;
% %    pos = [0.37 0.05 0.28 0.32];
% %    set(gcf,'Units','normalized');
% %    set(gcf,'Position',pos);
% %    set(gcf,'color','w'); 
% %       
% %    yyaxis left 
% %    xP = 1e3.*t; yP = vC;
% %    plot(xP,yP,'b','linewidth',2');
% %    grid on
% %    xlabel('time  t  [ms]');
% %    ylabel('v_C  [ V ]');
% %     
% %    ax = gca;
% %    ax.YAxis(1).Color = 'b';
% %    set(gca,'fontsize',FS);
% %    set(gca,'xLim',[0 1e3*max(t)]);
% %    
% %    yyaxis right
% %    xP = 1e3.*t; yP = 1e3.*iC;
% %    plot(xP,yP,'r','linewidth',2');
% %    ylabel('i_C  [ mA ]');
% %    ax = gca;
% %    ax.YAxis(2).Color = 'r';
% %    
% %    tm  = ' CAPACITOR: voltage phase lags current phase by \phi = -\pi / 2 rad'; 
% %    title(tm,'fontweight','normal','fontsize',12);
% % 
% % figure(3);   % INDUCTOR --------------------------------------------------
% %    FS = 14;
% %    pos = [0.67 0.05 0.28 0.32];
% %    set(gcf,'Units','normalized');
% %    set(gcf,'Position',pos);
% %    set(gcf,'color','w'); 
% %       
% %    yyaxis left 
% %    xP = 1e3.*t; yP = vL;
% %    plot(xP,yP,'b','linewidth',2');
% %    grid on
% %    xlabel('time  t  [ms]');
% %    ylabel('v_L  [ V ]');
% %     
% %    ax = gca;
% %    ax.YAxis(1).Color = 'b';
% %    set(gca,'fontsize',FS);
% %    set(gca,'xLim',[0 1e3*max(t)]);
% %    
% %    yyaxis right
% %    xP = 1e3.*t; yP = 1e3.*iL;
% %    plot(xP,yP,'r','linewidth',2');
% %    ylabel('i_L  [ mA ]');
% %    ax = gca;
% %    ax.YAxis(2).Color = 'r';
% %    
% %    tm  = ' INDUCTOR: voltage phase leads current phase by \phi = +\pi / 2 rad'; 
% %    title(tm,'fontweight','normal','fontsize',12);
% %    
% % figure(4)   % RESISTOR I/V curve ----------------------------------------
% %    FS = 14;
% %    pos = [0.07 0.45 0.28 0.32];
% %    set(gcf,'Units','normalized');
% %    set(gcf,'Position',pos);
% %    set(gcf,'color','w'); 
% %         
% %    xP = vR; yP = 1e3.* iR;
% %    plot(xP,yP,'b','linewidth',2');
% %    grid on
% %    xlabel('v_R [ V ]');
% %    ylabel('I_R  [ mA ]');
% %    set(gca,'fontsize',FS);
% %    %set(gca,'xLim',[0 1e3*max(t)]);
% %    tm = 'RESISTOR';
% %    title(tm,'fontweight','normal','fontsize',12);
% %    
% %  figure(5)   % CAPACITOR I/V curve ----------------------------------------
% %    FS = 14;
% %    pos = [0.37 0.45 0.28 0.32];
% %    set(gcf,'Units','normalized');
% %    set(gcf,'Position',pos);
% %    set(gcf,'color','w'); 
% %         
% %    xP = vC; yP = 1e3.* iC;
% %    plot(xP,yP,'b','linewidth',2');
% %    grid on
% %    xlabel('v_C [ V ]');
% %    ylabel('I_C  [ mA ]');
% %    set(gca,'fontsize',FS);
% %    %set(gca,'xLim',[0 1e3*max(t)]);
% %    tm = 'CAPACITOR';
% %    title(tm,'fontweight','normal','fontsize',12);  
% %    
% % figure(6)   % INDUCTOR I/V curve ----------------------------------------
% %    FS = 14;
% %    pos = [0.67 0.45 0.28 0.32];
% %    set(gcf,'Units','normalized');
% %    set(gcf,'Position',pos);
% %    set(gcf,'color','w'); 
% %         
% %    xP = vL; yP = 1e3.* iL;
% %    plot(xP,yP,'b','linewidth',2');
% %    grid on
% %    xlabel('v_L [ V ]');
% %    ylabel('I_L  [ mA ]');
% %    set(gca,'fontsize',FS);
% %    %set(gca,'xLim',[0 1e3*max(t)]);
% %    tm = 'INDUCTOR';
% %    title(tm,'fontweight','normal','fontsize',12);
% %    
% %    
% 
% 
% 
% % f1 = @(x,y) x^3 - 3*x*y^2 - 1;
% % f2 = @(x,y) 3*x^2*y - y^3;
% 
% % x = linspace(-1,1,1000);
% % y = linspace(-1,1,1000);
% % [X,Y] = meshgrid(x,y);
% % figure(1)
% % Z = (X.^3 - 3*X.*Y.^2 - 1).^2 + (3*X.^2.*Y - Y.^3).^2;
% % contour(X,Y,Z,50)
% % figure(2)
% % mesh(X,Y,Z)
% 
% 
% N = 24;
% h = 2/N;
% count = 100;
% figure;
% hold on;
% axis([-2,2,-2,2]);
% % red --- (0.8,0.1)
% % green --- (-0.5,0.8)
% % blue  ---- (-0.5,-0.8)
% for i = -N : N
%     x1 = i * h;
%     for j = -N : N
%         x2 = j * h;
%         x0 = [x1;x2];
%         if  x1 ~= 0 || x2 ~= 0
%             [x,iter] = newton(x0,10^-8,count);
%             if iter < count
%                 if x(1) > 0
%                     plot(x0(1),x0(2),'.','markersize',9,'color','red');
%                 elseif x(2) > 0
%                     plot(x0(1),x0(2), '.','markersize',9,'color', 'green');
%                 else 
%                     plot(x0(1),x0(2), '.','markersize',9,'color','blue');
%                 end
%             end
%             hold on;
%             drawnow;
%         end
%     end
% end
% 
% 
% % if x(2) < 0, form = '.r' .. at the end plot(x(1),x(2), '.r')
% 
%             
%           

% 
% N = 4;
% A = [1 0 0 N; 0 1 0 2*N; 0 0 1 3*N]
% B = A * A'
% 
% t = 0 : 0.01 : 10;
% A = 4;
% fc = 2;
% u = A * cos(2*pi*fc*t);
% kp = 20;
% kf = 20;
% a = 1;
% fm = 0.1;
% m = square(t*(0.2*pi));
% y = t .* (t < 5) + (10 - t).*(t >= 5);
% figure(1);
% plot(t,m);
% up = A * cos(2*pi*fc*t + kp * a * m);
% uf = A * cos(2*pi*fc*t + kf * a * y);
% subplot(2,1,1);
% plot(t,up);
% subplot(2,1,2);
% plot(t,uf);
% alpha = 0.905;
% JF = [2 + alpha/sqrt(alpha + 0.2^2); 2 + alpha/sqrt(alpha + 0.4^2);
%  2 + alpha/sqrt(alpha + 0.6^2);2 + alpha/sqrt(alpha +0.8^2);2 + alpha/sqrt(alpha + 1)^2]
% F = [ 2 * alpha + sqrt(alpha^2 + 0.2^2) - 1.55;
%   2 * alpha + sqrt(alpha^2 + 0.4^2) - 1.65;
%   2 * alpha + sqrt(alpha^2 + 0.6^2) - 1.8;
%   2 * alpha + sqrt(alpha^2 + 0.8^2) - 1.95;
%   2 * alpha + sqrt(alpha^2 + 1^2) - 2.1]
% A = JF' * JF
% b = -JF' * F

% A = [8,2,-2,1;2,9,2,9;-2,2,4,5;1,9,5,12]
% eig(A)

% format long
% A = [1 1; exp(1/2),1; exp(1),1]
% b = [-0.4,0.7,2.5]'
% A\b

% format rat
% A = [0.5,1.5,1.5;1,3,1;0,2,0]
% b = [10,5,6]'
% [U,V,p] = lu(A)
% x = A\b
% 
% x = -20:0.01:20;
% y = 2*cos(1/3*x);
% f1 = @(x) 2*cos(1/3*x);
% h1 = plot(x,y);
% hold on
% h2 = plot(x,x);
% x0 = 4;
% h3 = plot(x0,f1(x0),'o');
% for k = 1 : 10
%     x1 = f1(x0);
%     pause(0.2);
%     delete(h3);
%     h3 = plot(x1,f1(x1),'o');
%     x0 = x1;
% end
    
% p = 0.5;
% [X,Y] = meshgrid(-1:0.01:1);
% T = abs((X.^2 + Y.^2)).^p;
% Z = X./T;
% W = Y.*sin(1./X);
% U = X.*Y./(X.^2+Y.^2+eps);
% O = (X+Y).*sin(1./(X+eps)).*sin(1./(eps+Y));
% mesh(X,Y,O);
% 
% p = 0.85
% n = 6
% delta = (1-p)/n
% i = [2 3 4 4 5 6 1 6 1]
% j = [1 2 2 3 3 3 4 5 6]
% 
% G = sparse(i,j,1,n,n)
% full(G)
% c = full(sum(G))
% 1./c'
% D = spdiags(1./c',0,n,n)
% full(D)
% I = speye(6)
% x = (I - p*G*D)\(delta*e)
% A = p*G*D + delta
% x = (I-A)\e
% x = x/sum(x)
% e = ones(n,1);
% A = spdiags([e -2*e 3*e], -1:1, n, n)





% 
% 
% A = [1 2 3;4 5 6;7 8 9]
% b = [1 3 5]'
% det(A)
% x = A\b

% A = [2 -1 0 0 0 0 ; -1 2 -1 0 0 0; 0 -1 2 -1 0 0; 0 0 -1 2 -1 0; 0 0 0 -1 2 -1; 0 0 0 0 -1 2]
% x0 = 1/sqrt(1+4+9+16+25+36)*[1 2 3 4 5 6]'
% tol = 1e-10
% lambda0 = max(eig(A)); % groesster Eigenwert von A
% figure(5)
% hold on
% [lambda,x,k,v] = directVecIter(A,x0,tol,lambda0)
% plot(v);
% figure(6)
% hold on
% [lambda1,x1,k1,v1] = directVecIter(A,x0,1e-1,lambda0)
% [lambda2,x2,k2,v2] = inverseVecIter(A,x1,lambda1,tol,lambda0,k1,v1)
% plot(v2);
% figure(7)
% semilogy(abs(v2));

% A = [6,1,1,2;3,3,1,4;3,1,-4,2;1,2,3,6]
% v1 = [3,3]'+[3*sqrt(2),0]'
% Qv1 = eye(2) - 2*v1*v1'/(v1'*v1)
% Q1 = blkdiag(1,Qv1)
% AA = Q1*A*Q1
% [Q,R] = qr(AA)
% R*Q

% A = magic(10)
% s = Spektrum(A,1e-12)
% eig(A)

% V = zeros(3,22);
% for i = 4 : 1 : 25
%     
%     B = matrixB(i)
%     v = eig(B)
%     w = sort(v)
%     w1 = w(1:3)
%     V(1:3,i) = w1;
% end
% 
% plot(V') % the column of V' is the x -axis in which case is the n 
% % the row is the data at each index n (the n value on the x axis)
% xlim([4,26])
% 
% [V,D] = eig(matrixB(50));
% 
% eig_value = diag(D);
% [asort,ind] = sort(eig_value);
% Msort = V(:,ind);



% n = 36;
% temp = ones(n,1);
% A = diag(4*temp) - diag(1*temp(1:(n-1),1),-1)-diag(temp(1:(n-1),1),1) - diag(temp(1:(n-6),1),6) - diag(temp(1:(n-6),1),-6)
% b = [
%     20
%     10
%     10
%     10
%     10
%     20
%     10
%      0
%      0
%      0
%      0
%     10
%     10
%      0
%      0
%      0
%      0
%     10
%     10
%      0
%      0
%      0
%      0
%     10
%     10
%      0
%      0
%      0
%      0
%     10
%     20
%     10
%     10
%     10
%     10
%     20]
% x0 = zeros(n,1);
% r = b - A * x0;
% d = r;
% for k = 1 : 12
%     alpha = norm(r)^2 / ((A*d)'*d);
%     x_temp = x0 + alpha * d;
%     r_temp = r - alpha * A * d;
%     beta = norm(r_temp)^2 / norm(r)^2;
%     d = r_temp + beta * d;
%     r = r_temp;
%     x0 = x_temp;
% end
% x0
% x = A\b
% 
% norm(x-x0)/norm(x)
% 
% 
% dif = A*x0 - b
% 
% v = 10 * ones(36,1)
% result = A * v
% diff = A * v - b




% 

% 7.7 Beispiel 7.29
% A = [2.8452,1.5151,3.8814;1.3423,1.8106,2.8356;0.1438,-0.77 1.3433]
% A = [2,3,1;1,8,9;10,12,11];
% A = A' + A
% eig(A)
% 
% Qk = eye(3);
% for i = 1 : 100
%     B = A * Qk;
%     [Qk,Rk] = qr(B);
%     if mod(i,10) == 0
%         disp("Rk: ");
%         disp(Rk);
%         Ak = Qk' * A * Qk;
%         disp("Ak : ");
%         disp(Ak);
%         
%     end
% end
%     



% Beispiel 11.2

A = [1 1/2 0 -1; 1/2 -1 -1 0; 0 2 1 2; 2 0 -2 -1]
A1 = [1 1/2 0 -1; 0 1 2 -1; 0 0 -3 4; 0 0 3/2 -3/4]
A2 = [1 1/2 0 -1; 0 1 2 -1; 0 0 3 -4; 0 0 3/2 -3/4]
b = [0,0,-3,0]'
x = A\b
y = A1\b