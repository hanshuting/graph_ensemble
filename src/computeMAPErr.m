function err = computeMAPErr(F, G, Ut, Vt, YN, ixNode, ixEdge, edges)
% computeErrMAP  Compute Hamming error of the node MAP error
%
%   U, V, YN, YE, edges: same format as used in the Ising constructor.

    thetaN = F * Ut;
    thetaE = G * Vt;

    M = size(ixNode, 1);
    
    mistakes = 0;
    for m = 1:M
        ixNm = ixNode(m,1):ixNode(m,2);
        ixEm = ixEdge(m,1):ixEdge(m,2);
        
        [YhatNm, ~, eBelow] = ...
            solveQPBO(thetaN(:,ixNm), thetaE(:,ixEm), edges{m});
        
        mistakes  = mistakes + sum(vec(YhatNm' ~= YN(ixNm,:)));
    end
    
    % ixNode(end) = # of nodes.
    err = mistakes / (2*ixNode(end));
end

