function [M, cost, infeasible] = onematch(W)
% Compute perfect 1-matching of weighted adjacency matrix W.
%
%   [M, cost, infeasible] returns perfect matching as a sparse matrix M
%   and cost, or else sets infeasible to true.

    infeasible = false;

    [ri, rj, rw] = findUT(W);
    nNodes = max(max(ri), max(rj));
    edges  = blossom_double_mex(nNodes, ri, rj, rw);
    
    M = symmSparse(ri(edges), ...
                   rj(edges), ...
                   ones(sum(edges), 1), ...
                   nNodes);
                   
    cost = sum(W(logical(M(:))));
   

end