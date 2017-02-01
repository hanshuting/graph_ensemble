function [Aeq, beq] = makeOvercompleteLocalPolytope(W, L)
% makeOvercompleteLocalPolytope  Equality constraints for local polytope.

    N = size(W, 1);
    [iVec, jVec, wVec] = findUT(W);
    nEdges = length(iVec);
    
    % One constraint for each variable, and both sides of each edge per label    
    nConstr = N + 2*L*nEdges;
    nDim    = L*N + L*L*nEdges;
    
    % TODO: Sparsify later?
    Aeq(nConstr, nDim) = 0;    
    beq(nConstr, 1) = 0;
    beq(1:N)     = 1;      
    
    % Singleton marginalization constraints
    for n = 1:N
        % copy-on-write
        Aeq(n, ((n-1)*L+1):n*L) = 1;
    end
        
    cube(L, L, nEdges) = 0;
    rect(L, N) = 0;    
    
    cIdx = N + 1;
    for ii = 1:nEdges
        i = iVec(ii);
        j = jVec(ii);
                        
        % copy-on-write        
        % Convention: cube(k_i, k_j, ii)
        for k = 1:L
            % Marginalize out k_j; k_i is left.
%             cube(k,:,ii) = 1;
            % POSSIBLE MISTAKE
            cube(:,k,ii) = 1;
            
            rect(k,i)    = -1;
%             rect(k,j)    = -1;
            
            Aeq(cIdx,:) = packOvercomplete(rect, cube);
            cIdx = cIdx + 1;
            
            % clear                        
            cube(:,k,ii) = 0;
            rect(k,i)    = 0;
                        
            % Marginalize out k_i; k_j is left.
%             cube(:,k,ii) = 1;
            % POSSIBLE MISTAKE
            cube(k,:,ii) = 1;
            
            rect(k,j)    = -1;
%             rect(k,i)    = -1;
            
            Aeq(cIdx,:) = packOvercomplete(rect, cube);
            cIdx = cIdx + 1;
            
            % clear
            cube(k,:,ii) = 0;
            rect(k,j)    = 0;                        
        end       
    end
    
    assert(nConstr == cIdx - 1, 'Wrong number of constraints.');

end

