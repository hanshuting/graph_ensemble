function y = safelog(x);

if (ismember(0,x))
    y=zeros(size(x));
    y(find(x>0)) = log(x(find(x>0)));
    y(find(x==0)) = -inf;
else
    y = log(x);
end
%disp([num2str(x) num2str(y)]);