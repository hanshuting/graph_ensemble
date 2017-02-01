function g = ezgrad(f,x,order)

% calculate gradient by finite differences
if nargin < 3
    order = 2;
end

f0 = f(x);
e = sqrt(eps)*(1+norm(x(:),inf))*1000;
%e = 1e-7;

if order==1
    f2 = zeros(size(x));
    for i=1:numel(x);
        x2 = x;
        x2(i) = x(i) + e;
        f2(i) = f(x2);
    end
    g = (1/e)*(f2-f0);
elseif order==2
    f2 = zeros(size(x));
    f3 = zeros(size(x));
    for i=1:numel(x);
        x2 = x;
        x2(i) = x(i) + e;
        f2(i) = f(x2);
        x3 = x;
        x3(i) = x(i) - e;
        f3(i) = f(x3);
    end
    g = (1/e/2)*(f2-f3);
elseif order==-1
    e = sqrt(eps)*(1+norm(x(:),inf));
    %e = eps;
    g  = zeros(size(x));
    i  = sqrt(-1);
    for j=1:numel(x)
        x2 = x;
        x2(j) = x(j) + e*i;
        g(j) = imag(f(x2))/e;
    end
else
    error('order must be 1 or 2 or -1');
end