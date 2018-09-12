function M = csaAssignPermMat(C)
% csaAssignPerm  One-line wrapper for CSA

    N = size(C, 1);
    [spC, scale] = sparsifyAndRound(C - min(C(:)));
    edges = csaAssign(2*N, spC);
    perm = edges(2,:) - N;
    M = expandPerm(perm)';

end

