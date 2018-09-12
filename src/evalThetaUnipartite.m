function [hammingErr, meanAuc] = evalThetaUnipartite(theta, features, Ytrue)
%evalTheta  Evaluate normalized Hamming error of theta on features and Ytrue, cell Ytrue.
%
%   The normalized Hamming error is defined to be the fraction of edges in
%   the predicted matchings that are in the true matching.

M = length(Ytrue);
K = size(features, 2);

err = 0;
trials = 0;

sumAuc = 0;

for m = 1:M    
    % Construct specific A matrix for sample m  
    A = zeros(size(features{m,1}));
    for k = 1:K
        A = A + theta(k) .* features{m,k};
    end   

    [ai, aj, aw] = findUT(A);
    nNodes = max(max(ai), max(aj));
    edges  = blossom_double_mex(nNodes, ai, aj, -aw);

    medges = [ai(edges) aj(edges)];

    [yi, yj, ~] = findUT(Ytrue{m});
    yedges = [yi yj];

    assert(size(medges, 1) == size(yedges, 1), 'medges and yedges do not have the same number of rows');

    err    = err    + size(setxor(medges, yedges, 'rows'), 1);
    trials = trials + size(union( medges, yedges, 'rows'), 1);

    % meanAuc computation
    yVec = vecUT(Ytrue{m}, true);
    Avec = vecUT(A, true);

%    [x, y, t, meanAuc] = perfcurve(yVec, Avec, 1);
    sumAuc = sumAuc + fastAUC(yVec, Avec, 1);
end

% Normalized Hamming error
hammingErr = err / trials;
meanAuc    = sumAuc / M;
meanAuc
end

