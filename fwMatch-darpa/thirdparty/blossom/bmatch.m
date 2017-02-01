function [M, cost, infeasible] = bmatch(W, b, varargin)
% bmatch  Solve max-weight (perfect) variable b-matching
%   M = bmatch(W, b, perfect) returns a (sparse) matrix M indicating edges
%   chosen in the matching. W is an NxN (sparse) weighted graph; only upper
%   triangular matter. b is an Nx1 integer vector, or a scalar, indicating
%   number of matched needed for the node. perfect is a boolean controlling
%   whether equality or <= inequality is used for the b settings.
%
%   Code for imperfect matching uses the reduction in Section 1.5.1 of
%   Guido Schäfer, "Weighted Matchings in General Graphs", Master's Thesis, 
%   Universität des Saarlandes, 2000
    

    p = inputParser;
    p.KeepUnmatched = true; % Enable extraneous options.
    
    p.addRequired('W', @(WW) isnumeric(WW) && ...
        size(WW, 1) == size(WW, 2) && ...        
        all(W(:) >= 0));
    
    p.addRequired('b', @(bb) isnumeric(b) && ...
        length(bb) == size(W, 1));
    
    p.addParamValue('perfect', true);
    p.addParamValue('verbose', true);
    p.addParamValue('round', false);    
    
    p.parse(W, b, varargin{:});
    o = p.Results;
    
    nOrigNodes = size(W, 1);
                
    if o.verbose
        fprintf(1, 'bmatch with %d nodes, mean b = %g\n', size(W, 1), mean(b));
    end
        
    tic;
    [ri, rj, rw, em, periNodes, neg] = bondyMurty(-W, b, o.round);
    assert( (~neg && all(rw <= 0)) || (neg && all(rw >= 0)), 'Incoherent minimization problem');    
    bmTime = toc;
    
    if o.verbose
        fprintf(1, 'bondyMurty took %g secs\n', bmTime);
    end
        
    nNodes = max(max(ri), max(rj));    
    nEdges = length(ri);
    nBMEdges = nEdges;
            
    if ~o.perfect        
        tic;

        % Every vertex n in 1:nNodes is isomorphic to n + nNodes.
        shadowRi = ri + nNodes;
        shadowRj = rj + nNodes;
        
        if neg
            extraRi = find(~periNodes);
            
            % Translate up
            %rwNzSel = rw < 0;
            %rw(rwNzSel) = rw(rwNzSel) - min(rw) + 1;            
        else
            extraRi = find(periNodes);            
        end
                
        extraRj = extraRi + nNodes;
        nExtra  = length(extraRj);
        
        ri = horzcat(ri, shadowRi, extraRi);
        rj = horzcat(rj, shadowRj, extraRj);        
        rw = horzcat(rw, rw, zeros(1, nExtra));
        
        % Expand em from nNodes x nNodes to (2*nNodes) x (2*nNodes).
        % This requires copying the original em into the lower-right
        % nNodes x nNodes block.
        %em = blkdiag(em, em);
        
        nNodes = 2 * nNodes;
        nEdges = length(ri);
        
        perfTime = toc;
        if o.verbose
            fprintf(1, 'Imperfect reduction took %g secs\n', perfTime);
        end
    end
    
    % Bail out if we have an odd number of nodes
    infeasible = false;
    
    if mod(nNodes, 2) == 1
        M = sparse(nOrigNodes, nOrigNodes);
        if neg
            
            % But if we negated the graph, we're not really infeasible. We
            % should instead be full
            M = negFlip(W, M);
            cost = calcCost(W, M);            
        else
            warning('Infeasible!');            
            cost = [];
            infeasible = true;     
        end
        return;
    end
        
    tic;
    
    if o.round
         mVec = blossom_mex(nNodes, ri, rj, rw);        
    else
        mVec = blossom_double_mex(nNodes, ri, rj, rw);
    end
    blossomTime = toc;
    
    if o.verbose
        fprintf(1, 'blossom took %g secs\n', blossomTime);
    end
    
    % Select the expanded edges
    
    % In case we did the imperfect reduction, only select the edges in the
    % Bondy-Murty graph.
    origMVec = mVec(1:nBMEdges);
    
    % Peripheral edges have nonzerm em
    peripheral = em > 0;
    % Matched edges are both peripheral and belong to the original graph.
    matchEdges = em(peripheral & origMVec);
            
    M = symmSparse(real(matchEdges), imag(matchEdges), true(length(matchEdges), 1), nOrigNodes);
    
    if neg
        M = negFlip(W, M);
    end
    
    cost = calcCost(W, M); 
end

function M = negFlip(W, M)
    N = size(W, 1);
    W = triu(W) + triu(W)';

    % Relative complement        
    if issparse(W)
        M = logical(spones(W) - M);
    else
        pattern = ones(size(W)) - diag(ones(N, 1));
        M = logical(pattern - M);
    end
end

function c = calcCost(W, M)
    Wu = triu(W);
    Mu = triu(M);                

    c = sum(full(vec(Wu(Mu))));    
end
