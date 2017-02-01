function [Y, features] = makeDisconnectedCyclesBiased(theta, subD, M)
%makeBiasedDisconnectedCycles  Create permutations and dataset for disc. cycles
%
%   [Y, features] = makeBiasedDisconnectedCycles(theta, subD, M) makes a
%   (K*subD) x (K*subD) biadjacency matrix with a heavy diagonal term and
%   small full term each.

p = inputParser;

o = p.Results;

K = length(theta) / 2;
D = K * subD;

% Make random features.
for m = 1:M    
    for k = 1:K        
        fullK = k;
        diagK = K + k;        
        
        features{m,fullK} = zeros(D);
        features{m,diagK} = zeros(D);
        startIdx = subD * (k - 1) + 1;
        stopIdx  = startIdx + subD - 1;
        features{m,fullK}(startIdx:stopIdx,startIdx:stopIdx) = rand(subD);                        
        features{m,diagK}(startIdx:stopIdx,startIdx:stopIdx) = diag(rand(subD,1));                   
    end
end

% Now generate the true A matrix.
Y = zeros(M, D);
for m = 1:M
    Atrue = zeros(D);
    for k = 1:(2*K)
        Atrue = Atrue + theta(k) * features{m,k};
    end
    
    expA = exp(Atrue);
    expA(Atrue == 0) = 0;
    
    expA;
    
    Y(m,:) = sample_perms(expA, 1);
end

end

