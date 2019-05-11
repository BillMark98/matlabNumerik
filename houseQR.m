%% QR = houseQR (A) berechnet die QR-Zerlegung via Householder-Transformation
function QR = houseQR (A)

    [m,n]= size(A);
    
    for k= 1:n
        % Berechnung von v
        y= A(k:m,k);
        e= zeros(m-k+1,1);
        e(1)= 1;
        alpha= norm(y);
        if (y(1)<0)
            % Ausloeschung vermeiden
            alpha= -alpha; 
        end    
        v= y + alpha*e;
        
        disp 'alpha = '
        disp (alpha)
        disp '[y, v]= '
        disp ([y, v])
        
        % Anwenden von Q_v:
        % Q_v y = -alpha e. Wir speichern die Nullen nicht, da sie gleich
        % mit dem skalierten v ueberschrieben werden.
        A(k,k)= -alpha;
        % Q_v in Spalten j>k anwenden
        beta= 2/(v'*v);
        for j= k+1:n
            z= A(k:m,j);
            A(k:m,j)= z - beta*(v'*z)*v;
        end
        % skaliertes v unterhalb der Diagonalen abspeichern
        v= v/v(1);
        A(k+1:m,k)= v(2:end);
    end
    
    QR= A;
    
end
        
