
% uebung 05

% f = @(x) x.^2 + 2*x + x.^3 % 2/5*(5/2*x^3 - 3/2*x) + 2/3*(3/2*x^2 - 1/2) + 3/5*x+2*x + 1/3
% % = 2/5 * phi3 + 2/3 * phi2 + 2.6 * phi1 + 1/3
% u = legendreL2Approximation(20,f);
% 
% v = monomL2Approximation(5,f)

% u = [1,2,1]'
% x = [0,1]'
% % y = 1/2 + 2x + 3/2*x^2
% y = legrendrePolynom(u,x)
% 
% u1 = [1,1,2,2]'
% x1 = [0,0.5,1]'
% % y1 = 5x^3 + 3x^2 - 2x
% y1 = legrendrePolynom(u1,x1)


L = [10,50,100];
f1 = @(x) abs(x)
for i = 1 : 3
    figure(i);
% 
%     subplot(2,1,1);    
%     hold on;
%     pnM = monomL2Approximation(L(i),f1);
%     pnL = legendreL2Approximation(L(i),f1);
%     fplot(f1,[-1,1],'color','blue');
%     fplot(@(x)(monomPolynom(pnM,x)),[-1,1],'color','magenta');
%     legend('f(x)','monom');
%     
%     subplot(2,1,2);
%     hold on;
%     fplot(f1,[-1,1],'color','blue');
%     fplot(@(x)(legrendrePolynom(pnL,x)),[-1,1],'color','cyan');
%     legend('f(x)','legendre');
%     legend('f(x)','monom','legendre');


         
    hold on;
    pnM = monomL2Approximation(L(i),f1);
    pnL = legendreL2Approximation(L(i),f1);
    fplot(f1,[-1,1],'color','blue');
    fplot(@(x)(monomPolynom(pnM,x)),[-1,1],'color','magenta');
    legend('f(x)','legendre');
    fplot(@(x)(legrendrePolynom(pnL,x)),[-1,1],'color','cyan');
    legend('f(x)','monom','legendre');
end