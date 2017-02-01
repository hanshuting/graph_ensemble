function [x,o] = gradDescent(fun,x,opts)

for t = 1:opts.maxIters
    [val, grad] = fun(x);
    x = x - opts.stepSize(t)*grad;    
end
[o, grad] = fun(x);
