% x= lr_solve(LR,p,b) berechnet zu gegegebener LR-Zerlegung A=L*R und 
% Zeilenpermutation p die Loesung von A*x = b durch
% Vorwaerts-/Rueckwaertseinsetzen.
function x= lr_solve( LR, p, b)

    n= length(b);
    y= zeros(size(b));
    x= zeros(size(b));
    
    % permutiere b
    b= b(p);
    
    % loese L*y = b
    for i=1:n
        y(i)= b(i) - LR(i,1:i-1)*y(1:i-1);
    end

    % loese R*x = y
    for i=n:-1:1
        x(i)= (y(i) - LR(i,i+1:n)*x(i+1:n))/LR(i,i);
    end

end

