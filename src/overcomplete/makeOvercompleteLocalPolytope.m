function [Aeq, beq] = makeOvercompleteLocalPolytope(N, L, pairs)
% makeOvercompleteLocalPolytope  Equality constraints for local polytope.

%     N = size(W, 1);
%     [iVec, jVec, ~] = findUT(W);
    nEdges = size(pairs, 1);
    
    % One constraint for each variable, and both sides of each edge per label    
    nConstr = N + 2*L*nEdges;
    nDim    = L*N + L*L*nEdges;
    
    nnz     = N*L + 2*nEdges*L;
    
    Ai = zeros(nnz, 1);
    Aj = zeros(nnz, 1);
    Aw = zeros(nnz, 1);
        
    % TODO: Sparsify later?
%     Aeq(nConstr, nDim) = 0;    
    beq(nConstr, 1) = 0;
    beq(1:N)     = 1;      
    
    % Singleton marginalization constraints
    for n = 1:N
        js = ((n-1)*L+1):n*L;
        Ai(js) = n * ones(L, 1);
        Aj(js) = js;
        Aw(js) = 1;
    end
        
    nodeLookup = reshape(1:L*N, L, N);
    edgeLookup = reshape(1:L*L*nEdges, L, L, nEdges);
    
    % Shift all indices up to begin at the edge indices; the first L*N
    % indices are for nodes.
    edgeLookup = edgeLookup + L*N;
    
    Aidx = N*L + 1;
    cIdx = N + 1;
    for ii = 1:nEdges
        i = pairs(ii,1);
        j = pairs(ii,2);
                        
        % copy-on-write        
        % Convention: cube(k_i, k_j, ii)
        for k = 1:L
            % Marginalize out k_j; k_i is left.
            
            Ai(Aidx) = cIdx;
            Aj(Aidx) = nodeLookup(k,i);
            Aw(Aidx) = -1;
            
            Aidx = Aidx + 1;
            edgeAidxs = Aidx:Aidx+L-1;
            Ai(edgeAidxs) = cIdx;
            Aj(edgeAidxs) = edgeLookup(:,k,ii);
            Aw(edgeAidxs) = 1;
            
            Aidx = edgeAidxs(end) + 1;
            cIdx = cIdx + 1;
                        
            % Marginalize out k_i; k_j is left.
            Ai(Aidx) = cIdx;
            Aj(Aidx) = nodeLookup(k,j);
            Aw(Aidx) = -1;
            
            Aidx = Aidx + 1;
            edgeAidxs = Aidx:Aidx+L-1;
            Ai(edgeAidxs) = cIdx;
            Aj(edgeAidxs) = edgeLookup(k,:,ii);
            Aw(edgeAidxs) = 1;
            
            Aidx = edgeAidxs(end) + 1;
            cIdx = cIdx + 1;
        end       
    end
    
    assert(nConstr == cIdx - 1, 'Wrong number of constraints.');

    Aeq = sparse(Ai, Aj, Aw, nConstr, nDim);
    
end

