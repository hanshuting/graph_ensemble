function [Y, features] = makeDisconnectedCycles(theta, subD, M, varargin)
%makeDisconnectedCycles  Create permutations and dataset for disc. cycles
%
%   [Y, features] = makeDisconnectedCycles(theta, subD, M) makes a
%   (K*subD) x (K*subD) biadjacency matrix 

p = inputParser;

p.addRequired('theta', @isnumeric);
p.addRequired('subD',  @isscalar);
p.addRequired('M',     @isscalar);
p.addParamValue('subGenerator', @rand);

p.parse(theta, subD, M, varargin{:});

o = p.Results;

K = length(theta);
D = K * subD;

% Make random features.
for m = 1:M    
    for k = 1:K
        features{m,k} = zeros(D);
        startIdx = subD * (k - 1) + 1;
        stopIdx  = startIdx + subD - 1;
            
        features{m,k}(startIdx:stopIdx,startIdx:stopIdx) = o.subGenerator(subD);
    end
end

% Now generate the true A matrix.
Y = zeros(M, D);
for m = 1:M
    Atrue = zeros(D);
    for k = 1:K
        Atrue = Atrue + theta(k) * features{m,k};
    end
    
    expA = exp(Atrue);
    expA(Atrue == 0) = 0;
    
    expA
    
    Y(m,:) = sample_perms(expA, 1);
end

end

