function [ T ] = randTree( nNodes, excessDensity )
% T = randTree(nNodes, nLoops)
%
% Create a random uniform spanning tree T of nNodes nodes with Wilson's
% algorithm. Optionally add nLoops random edges.
%
% Simplified for unweighted complete graph.
%
% T is a sparse nNodes x nNodes adjacency matrix.
%
% References
% [1] David Bruce Wilson, "Generating Random Spanning Trees More Quickly
%     than the Cover Time" STOC '96 doi>10.1145/237814.237880

    if nargin == 1
        nLoops = 0;
    end

    r = randsample(nNodes, 1);
    parents = zeros(nNodes, 1);
    inTree = false(nNodes, 1);
    inTree(r) = true;
    for i = 1:nNodes
        u = i;
        while ~inTree(u)
            pop = 1:nNodes;
            pop(u) = [];
            ri = randi(nNodes - 1);
            parents(u) = pop(ri);
            u = parents(u);
        end
        u = i;
        while ~inTree(u)
            inTree(u) = true;
            u = parents(u);
        end
    end
    
    % the root does not have a parent    
    parentsIvec = find(parents);
    parentsJvec = parents(parentsIvec);
    wvec        = true(nNodes - 1, 1);
    
    T = sparse(parentsIvec, parentsJvec, wvec, nNodes, nNodes);    
    Wexcess = sprand(~T) < excessDensity;

    T = T + Wexcess;
    T = T | T';
end
