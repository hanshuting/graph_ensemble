function [errs, denoms, aucs] = evalThetaUnipartiteSubsample(...
    theta, features, YTrue, nEdges, nTrials)
% evalThetaUnipartiteSubsample  Normalized Hamming error on a subsampled matrix
%
% NOTE: This version only tests one sample set. (trials refers to number of
%       random subsamples of a single sample of a matching.)
%
%       This means Ytrue is a sparse matrix, NOT a cell.

K = length(features);

A = zeros(size(features{1}, 1));
for k = 1:K
    A = A + theta(k) .* features{k};
end   

errs   = zeros(nTrials, 1);
denoms = zeros(nTrials, 1);
aucs   = zeros(nTrials, 1);

for t = 1:nTrials
    % Construct specific A matrix for sample m.
    % No references to the full A or YTrue beyond this point.
    [Asub, YTrueSub] = sampleSubmatching(A, YTrue, nEdges);    
    
    [ai, aj, aw] = findUT(Asub);
    nNodes = max(max(ai), max(aj));
    edges  = blossom_double_mex(nNodes, ai, aj, -aw);

    medges = [ai(edges) aj(edges)];

    [yi, yj, ~] = findUT(YTrueSub);
    yedges = [yi yj];

    assert(size(medges, 1) == size(yedges, 1), 'medges and yedges do not have the same number of rows');

    % KT 2/1 -- CHANGED to simply be
    % err:   number of edges in ytrue not predicted by the matching
    % denom: number of edges in ytrue
    %
    % The old definition double-counted and double-divided, which I think
    % is equivalent (I'm not positive though), but makes it harder to
    % directly compare to the normalization used in
    % randomMatchingHammingErr.m.
    errs(t)   = size(setdiff(yedges, medges, 'rows'), 1);
    denoms(t) = size(yedges, 1);

    % meanAuc computation
    yVec = vecUT(YTrueSub, true);
    Avec = vecUT(Asub, true);

    aucs(t) = fastAUC(yVec, Avec, 1);
end

end

