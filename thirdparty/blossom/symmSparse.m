function W = symmSparse(i, j, w, nNodes)
% symmSparse  Create symmetric sparse matrix from triangular factors.
%   [ W ] = symmSparse( iVec, jVec, wVec )
%
%   iVec, jVec, wVec - triangular part. Should be same dimensions.
%   W                - sparse matrix

    % Assume that all i, j, w share the same layout!
    d = find(size(i) > 1, 1);
    % Special case: singleton
    if isempty(d)
        W = sparse(i, j, w, nNodes, nNodes);
    else    
        W = sparse(cat(d, i, j), cat(d, j, i), cat(d, w, w), nNodes, nNodes);
    end

end

