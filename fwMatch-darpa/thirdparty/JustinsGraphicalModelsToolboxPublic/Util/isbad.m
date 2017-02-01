function o = isbad(x)
% tells you if x is bad (inf, nan)

xsum = sum(x(:));

o = isinf(xsum) || isnan(xsum);