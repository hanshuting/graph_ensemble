function [M, C] = blossom(A)
% blossom  Find min-cost perfect matching.
%    [M, C] = blossom(A) computes min-cost perfect matching on an undirected
%    graph weighted adjacency matrix A. Only the upper-triangle is kept.
%
%    M contains the returned matching, represented as a sparse adjacency matrix.
%    C is the cost.

    tic;
    nNodes = size(A, 1);
    [iVec, jVec, wVec] = findUT(A);
    eVec = blossom_double_mex(nNodes, iVec, jVec, wVec);

    M = sparse(iVec, jVec, eVec, nNodes, nNodes);
    C = sum(wVec(eVec));

   % fprintf('Total calltime: %g\n', toc);


end
