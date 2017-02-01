function [rhoNMat, rhoEMat] = makeRhoMatOvercomplete(edgeRhos, edges, N, L)
% makeRhoVecsOvercomplete  Make matricized rhos for a single sample.
%
%   rho = makeRhoOvercomplete(edgeRhos, edges, N, L)
%
%   edgeRows : E-vector of a reweighting parameter per edge.
%   N, L     : Number of nodes and labels
%
%   rhoNMat, rhoEMat : row vectors of rhos, to be s
%
%   rho      : Overcomplete vector of reweighting parameters

    assert(iscolumn(edgeRhos), 'edgeRhos must be a column vector.');

    E = length(edgeRhos);
    
    % Add the edge rhos to nodes on both sides of the edge.
    nodeRhos = zeros(N, 1);        
    for e = 1:E
        ix = edges(:,e);
        nodeRhos(ix) = nodeRhos(ix) + edgeRhos(e);
    end
    
    % Vectorization-foo to conform to the correct dimensions
    rhoNMat = repmat(nodeRhos, 1, L);                 % LxN
    rhoEMat = repmat(edgeRhos, 1, L^2);               % L^2 x E
end

