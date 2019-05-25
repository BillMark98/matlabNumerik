function [y] = aitkenNevilleIP(t,ft,x)
    n = length(t);
%     t = makeRow(t);
%     ft = makeRow(t);
    x = makeRow(x);
    m = length(x);
    % don't use recursion version because it's exponential complexity
    
%     if n == 1
%         y = ft(1) * ones(1,m);
%         return;
%     end
%     % remember have to write the whole formula in one line if write 
%     % y = (x - t(1))/(t(n) - t(1)) .* aitkenNevilleIP(t(2:end),ft(2:end),x)
%     % + (t(n) - x)/(t(n) - t(1)) .*
%     % aitkenNevilleIP(t(1:(end-1)),ft(1:(end-1)),x); will only evalute the
%     % first summand 
%     y = (x - t(1))/(t(n) - t(1)) .* aitkenNevilleIP(t(2:end),ft(2:end),x) + (t(n) - x)/(t(n) - t(1)) .* aitkenNevilleIP(t(1:(end-1)),ft(1:(end-1)),x);

      % use iteration instead
      for index = 1 : m
          A = zeros(n,n);
          A(:,1) = ft;
          for i = 2 : n
              for j = 2 : i
                  A(i,j) = A(i,j-1) + (x(index) - t(i))/(t(i) - t(i-j+1)) * (A(i,j-1) - A(i-1,j-1));
              end
          end
          y(index) = A(n,n);
      end
end
% convert a given vector to a row vector
function t = makeRow(t)
    n = length(t);
    [a,b] = size(t);
    if a == n
        t = t'
    end
end

