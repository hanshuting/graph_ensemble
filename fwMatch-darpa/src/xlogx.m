function y = xlogx(x)
if(x == 0)
    y =  0;
elseif(x == 1)
    y=  0;
else
    y = log(x^x);
end