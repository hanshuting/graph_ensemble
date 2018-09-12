function [ri, rj] = matUT(v)
% matUT  Matricize an upper triangular weight vector.
%
%   [ri, rj] = matUT(v) returns row and column indices in ri and rj
%   corresponding to the edges of a graph whose upper triangular weight
%   vector is given by v.

nEdges = length(v);
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


end