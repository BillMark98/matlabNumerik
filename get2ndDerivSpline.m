function m = get2ndDerivSpline(f,h)
f = MyTranspose(f);
n = length(f);

A = 4 * eye(n-2) + diag(ones(n-3,1),1) + diag(ones(n-3,1),-1);
v = f( 1 : end-2,1) - 2 * f(2 : end - 1,1) + f(3:end,1);
v = 6/(h^2) * v;
m = A\v;
m = [0;m;0];
end

% transpose y to column form if necessary
function y = MyTranspose(f)
[m,n] = size(f)
if(n > m && m == 1)
    y = f';
else
    y = f;
end

end