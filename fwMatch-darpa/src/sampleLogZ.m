function [logZ] = exactLogZRyser(features, theta)
    
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
    
    % And replace with -inf. Yay hacking math.
    logA(imag(logA) == 1) = -inf;
       
    A = exp(logA);    
    Z = ryser_small(A);
    assert(isfinite(Z));
    logZ = log(Z);
end

