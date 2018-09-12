function x2 = cummin(x1)
x2 = x1;
for i=2:length(x2)
    if x2(i-1) < x2(i)
        x2(i) = x2(i-1);
    end
end