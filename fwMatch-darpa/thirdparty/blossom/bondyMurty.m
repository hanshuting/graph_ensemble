function [ri, rj, rw, em, periNodes, neg] = bondyMurty(W, b, roundMe)
% bondyMurty  Find variable b-matching by reducing to unipartite matching.
%
%   [ri, rj, rw] = bondyMurty(W, b1) produces a unipartite matching problem
%   represented as a weighted graph with edge list (ri, rj) and weights rw,
%   some of which may be zero, using the procedure in Bondy & Murty page
%   431. See also Gerards 1995, Tutte (1954b).
%
%   E is a sparse complex matrix with same shape as the [ri, rj] graph.
%   E(i,j) is nonzero if the edge (i,j) in the expanded graph corresponds
%   to a node in the original W. The real part is the number of the lower
%   node in the original graph.
%
%   W  - (sparse) NxN matrix; only upper triangular matters.
%   b1 - Nx1 vector of matchings.
%
%   ri, rj - Edge list of extended graph. (ri(n), rj(n)) defines an edge.
%            Denote by nEdges = length(ri).
%   rw     - Weights of graph; can be zero.
%   periNodes - nNodes logical vector of peripheral nodes.
%   em     - nEdges complex vector. A peripheral edge in the expanded graph
%            has a complex a + bi entry mapping to edge (a, b) in the
%            original graph. Non-peripheral edges map to zero.
%
% Example:
% n=6;Q=rand(n,n);for i=1:n;Q(i,1:i)=0;end;[P]=edmondsreduction(Q,2);
% Example in Bondy&Murty
% Q=zeros(6,6);Q(1,2)=1;Q(1,4)=1;Q(1,6)=1;Q(2,3)=1;Q(3,4)=1;Q(3,5)=1;Q(4,5)=1;Q(4,6)=1.1;Q(5,6)=1;
% Solution is P(1,2);P(1,6);P(2,3);P(3,5);P(4,5);P(4,6);
%
%
%     0     1     0     0     0     1
%     1     0     1     0     0     0
%     0     1     0     0     1     0
%     0     0     0     0     1     1
%     0     0     1     1     0     0
%     1     0     0     1     0     0
%
% Qo=rand(n,n);for i=1:n;Qo(i,1:i)=0;end;P2=edmonds2(Qo,2);P=edmondsreduction(Qo,2);[sum(sum(P.*Qo)) sum(sum(P2.*Qo))]
%
% Another test, reverse weights should give same answer
% 
% n=6; b=2;
% Q=rand(n,n);for i=1:n;Q(i,1:i)=0;end;tic;P=edmondsreduction(Q,b); toc;
% nQ=max(max(Q))-Q;for i=1:n;nQ(i,1:i)=0;end;tic;nP=edmondsreduction(nQ,n-b-1);toc;
% nnP=1-nP;for i=1:n;nnP(i,1:i)=0;end;nnP=max(nnP,nnP');imagesc(nnP-P);
%
% If b<n/2, we do the reverse-weight solution since it's faster!

    if nargin < 3
        roundMe = true;
    end

    N = size(W, 1);      
    neg = false;
    
    % Solving a maximization problem over nonnegative weights. This is
    % equivalent to minimizing the negative weights.
    
    if issparse(W)    
        [iVec, jVec, wVec] = findUT(W);
    else
        % Keep all the weights. We may have zeros; keep them anyways                
        
        % Adding an irrational will do the trick to densify it.
        % No it doesn't... stupid numerical shit
%         er = eps(range(W(:)));
%         
%         nEls = N/2 * (N - 1);
%         iVec(nEls) = 0;
%         jVec(nEls) = 0;
%         wVec(nEls) = 0;
%         ind = 1;
%         for i = 1:N
%             for j = (i+1):N
%                 iVec(ind) = i;
%                 jVec(ind) = j;
%                 wVec(ind) = W(i,j);
%                 ind = ind + 1;
%             end
%         end
        
        [iVec, jVec, wVec] = findUT(W + eps);
        wVec = wVec - eps;
    end
        
    assert(all(wVec <= 0) || all(wVec >= 0), 'Edge weights must be monotonic');
    
    % BlossomV needs integers.
    %newRange = 2e8; % Magic number    
%         
%     % Scale everything, but don't mess with zero weights!
%     %sel = wVec ~= 0;
%     %assert(range(wVec(sel)) > 0);
%     
%     shiftScaleWVec = newRange * (wVec - max(wVec)) / (range(wVec) + 1) - 1;        
%     wVec = int32(shiftScaleWVec);    
%     assert(max(wVec) == -1);
    %assert(min(wVec) == -newRange - 1); 
    
%     fprintf(1, 'wVec diagnostic: min = %g, max = %g, range = %g\n', ...
%             min(wVec), max(wVec), range(wVec));            

    % Somehow just the 0.1 works better... much better.
    
    if roundMe
        if range(wVec) > 1
            wVec = wVec / range(wVec);        
        end

        swVec = 0.1 * double(intmax) * wVec; 
        iwVec  = int32(floor(swVec));    

        % Precision calculations
        uw  = length(unique(wVec));
        usw = length(unique(swVec));
        uiw = length(unique(iwVec));
        fprintf(1, 'Precision diagnostic: unique originals = %d, scaled = %d, integral = %d\n', ...
                uw, usw, uiw);

        % Precision calculations for the non-negligible
        thresh = 0.5;
        idxs   = abs(wVec) > thresh;
        uwt  = length(unique(wVec(idxs)));
        uswt = length(unique(swVec(idxs)));
        uiwt = length(unique(iwVec(idxs)));
        fprintf(1, 'Precision diagnostic at threshold %g: unique originals = %d, scaled = %d, integral = %d\n', ...
                thresh, uwt, uswt, uiwt);    
    else
        iwVec = wVec;
    end    
        
    % The number of times a node shows up in either iVec or jVec is its
    % degree.
    deg = histc(iVec, 1:N) + histc(jVec, 1:N);
    
    nCore = deg - b;      % nCore is |X_v| in Bondy and Murty.    
    if any(nCore < 0)
        error('bondyMurty:b', 'b must be less than deg for each node');
    end
    
    % Each vertex u in the original graph expands to deg(u) - b(u) core
    % vertices. If instead b(u) < deg(u) - b(u), we find a minimum-weight
    % matching of deg(u) - b(u) nodes (minimum weight in the original
    % problem, so we negate the weights again).    
    %
    % THIS DOES NOT WORK FOR IMPERFECT MATCHING! You need to derive a >= 1
    % matching, or a reduction from >= to perfect. Is this even possible?
    if (sum(b) < sum(nCore))        
        % Minimize the original problem (negate the weights once again)
        iwVec = -iwVec;
        nCore = b;
        b = deg - b;
        neg = true;
%         disp('Solving negated problem');
    end    
    
    % The X_v and Y_v sets contribute nodes.
    nVerts = sum(nCore) + sum(deg);
    % The edges in the bipartitions H_v = X_v \cup Y_v contribute edges. And
    % the number of edges in a complete bipartite graph is |X_v||Y_v|.
    nEdges = length(iVec) + sum(nCore .* deg);

    % Map each node to complete bipartite graph. Use the notation of Bondy and
    % Murty. This is just assigning unique node numbers.
    sG(N) = struct('core', [], 'peri', []);

    % Begin our extra vertices one after the end.
    ne = 1;
    ri = zeros(1, nEdges);
    rj = zeros(1, nEdges);
    
    if roundMe
        rw = zeros(1, nEdges, 'int32');
    else
        rw = zeros(1, nEdges);
    end
    em(nEdges) = 0;

    periNodes = false(1, nVerts);
    
    ind = 1;
    for n = 1:N
        % Make the bipartite nodes. This is really just allocating ids.
        sG(n).core = ind:(ind + nCore(n) - 1);
        sG(n).peri = (ind + nCore(n)):(ind + nCore(n) + deg(n) - 1);
        periNodes(sG(n).peri) = true;

        % Add edges to make the bipartite graph. They have weight zero, so they
        % don't affect the solution weight.
        %
        % MAY 6 -- VECTORIZED THIS SHIT!
        nAddEdges = length(sG(n).core) * length(sG(n).peri);        
        if nAddEdges > 0        
            [x, y] = meshgrid(sG(n).core, sG(n).peri);        
            nextNE = ne + nAddEdges - 1;
            ri(ne:nextNE) = x(:);
            rj(ne:nextNE) = y(:);
            rw(ne:nextNE) = 0;

            ne = nextNE + 1;
        end

        ind = ind + nCore(n) + deg(n);
    end
    assert(ind - 1 == nVerts);
    
    % Map each edge of the original graph onto a periphery edge in the
    % expanded graph. Remember its origin in em (the edge direction is
    % stored as a complex).    
    
    for e = 1:length(iVec)        
        i = iVec(e);
        j = jVec(e);                        
        
        % Mark a periphery vertex as "used" by setting it to zero. The
        % .peri are row vectors.
        [~, iPeri, vi] = find(sG(i).peri, 1);
        [~, jPeri, vj] = find(sG(j).peri, 1);

        sG(i).peri(iPeri) = 0;
        sG(j).peri(jPeri) = 0;

        ri(ne) = vi;
        rj(ne) = vj;
        rw(ne) = iwVec(e);
        em(ne) = complex(i, j);
                
        ne = ne + 1;
    end
    
    assert(ne - 1 == nEdges);    
end
