function x2 = polyexpand(x,order,x0)

% x2 = polyexpand(x,j)
% computes a polynomial basis expansion of order j
% j can be a vector in which case it computes several basis expansions and
% concatenates them.



if nargin < 3
    x0 = x;
end

if order==1
    x2 = x;
    return;
end

if length(order)>1
    x2 = [];
    for j=order(:)'
        x2 = [x2; polyexpand(x,j)];
    end
    return
end

z = [];
for i=1:size(x0,1)
    z = [z; x.*repmat(x0(i,:),size(x,1),1)];
end

x2 = polyexpand(z,order-1,x0);
