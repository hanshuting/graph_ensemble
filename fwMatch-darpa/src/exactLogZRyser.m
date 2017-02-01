function [logZ, dLogZ] = exactLogZRyser(features, theta)

    K = length(features);
    D = size(features{1}, 1);

    logA = zeros(D);
    % Sentinel
    
    if any(isinf(theta))
        warning('exactLogZRyser:theta', 'theta is not a real vector. Make you know what you''re doing!');
        theta(theta == -inf) = complex(0,1);
    end
    
    for k = 1:K
        logA = logA + theta(k) * features{k};
    end
    
    % And replace imaginary with -inf. Yay hacking math.
    logA(imag(logA) ~= 0) = -inf;
       
    A = exp(logA);
    assert(all(vec(isfinite(A))) && all(vec(isreal(A))));

    Z = ryser_small(A);
    assert(isfinite(Z));
    logZ = log(Z);

    dLogZ = zeros(D);
    for i = 1:D
        for j = 1:D
            expMe = A;
            expMe(i,j) = 0;
            ZMe = ryser_small(expMe);
            dLogZ(i,j) = 1 - ZMe / Z;
        end
    end        
   
    dLogZ = vec(dLogZ);
    
end

