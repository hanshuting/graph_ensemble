function rho = makeRhoOvercomplete(edgeRhos, edges, N, L)
% makeRhoVec  Make vectorized rhos
%
%   rho = makeRhoOvercomplete(edgeRhos, edges, N, L)
%
%   edgeRows : E-vector of a reweighting parameter per edge.
%   N, L     : Number of nodes and labels
%
%   rho      : Overcomplete vector of reweighting parameters

    assert(iscolumn(edgeRhos), 'edgeRhos must be a column vector.');

    E = length(edgeRhos);
    rho = zeros(N*L + E*L^2, 1);
    
    % Add the edge rhos to nodes on both sides of the edge.
    nodeRhos = zeros(N, 1);        
    for e = 1:E
        ix = edges(:,e);
        nodeRhos(ix) = nodeRhos(ix) + edgeRhos(ix);
    end
    
    % Vectorization-foo to conform to the correct dimensions
    nodeRhoMat = repmat(nodeRhos', L, 1);                 % LxN
    edgeRhoMat = repmat(edgeRhos', L^2, 1);               % L^2 x E

    rho = packOvercomplete(nodeRhoMat, edgeRhoMat);
end

