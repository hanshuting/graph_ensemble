function z = safeplus(x,y);

z = x+y;
z(find(isnan(z))) = 0;