function [m, cost] = vecUniMatch(w)
% vecUniMatch   Vectorized unipartite 1-matching.
%
%   [m, cost] = vecUniMatch(w) returns min-cost matching matrix m and its
%   cost given weight matrix m, where all matrices are given
%   as vectors of upper triangular entries.
%
% See also vecUT.

nEdges = length(w);
% Solve n*(n-1)/2 = c for c
nNodes = 0.5 * (sqrt(8*nEdges + 1) + 1);

% Fill in the row and column labels.
ridx = 1;
ri = zeros(nEdges, 1);
rj = zeros(nEdges, 1);

for col = 2:nNodes
    lastRow = col - 1;
    rrange  = ridx:(ridx + lastRow - 1);
    
    ri(rrange) = 1:lastRow;
    rj(rrange) = col;
    
    ridx = rrange(end) + 1;
end

m = blossom_double_mex(nNodes, ri, rj, w);

if nargout >= 2
    cost = sum(w(m));
end
end