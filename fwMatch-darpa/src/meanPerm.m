function mu = meanPerm(X)
% meanPerm  Compute empirical marginals of permutations in X

    if iscell(X)
        M = length(X);
        D = size(X{1}, 1);
        mu = zeros(D);
        for m = 1:M
            mu = mu + X{m};
        end
        mu = mu / M;
    else
        [M, D] = size(X);
        mu = zeros(D);
        for m = 1:M
            xp = expandPerm(X(m,:), @zeros);
            assert(all(sum(xp, 1) == 1) && all(sum(xp, 2) == 1));            
            mu = mu + expandPerm(X(m,:), @zeros);
        end
        mu = mu / M;
    end

end

