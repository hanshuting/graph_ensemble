function g = complexStepGrad(f, x, h)

    N = length(x);
    hVec = zeros(N, 1);
    g = zeros(N, 1);
    
    fx = f(x);
    
    for n = 1:N
        hVec(n) = h*j;
        g(n) = imag(f(x + hVec) / h);
        hVec(n) = 0;
    end


end

