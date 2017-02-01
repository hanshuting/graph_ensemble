function [ixNode, ixEdge] = computeUVIx(Ns, nEdges)
% computeUVIx  Compute begin and end row indices for a vectorized MRF.
%
%   ixNode(m,1) contains first row and ixNode(m,2) contains last row of U,
%   TN, and YN corresponding to sample m.

    M = length(Ns);
    assert(M == length(nEdges), 'Ns and nEdges must have same length.');
    
    ixNode = zeros(M, 2);
    ixEdge = zeros(M, 2);
    
    NRow = 1;
    ERow = 1;
    
    for m = 1:M
        ixNode(m,:) = [NRow, NRow + Ns(m) - 1];
        NRow = ixNode(m,2) + 1;
        
        ixEdge(m,:) = [ERow, ERow + nEdges(m) - 1];
        ERow = ixEdge(m,2) + 1;
    end

end

