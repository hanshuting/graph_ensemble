function [xNode, xEdge, ee] = solveMAP(gradNodes, gradEdges, Aeq, beq, L)

    [L, N] = size(gradNodes);
    [Lsq, nEdges] = size(gradEdges);
    
    assert(L*L == Lsq, 'gradNodes, gradEdges incompatible sizes');
    model.A = Aeq;
    model.rhs = beq;
    model.obj = packOvercomplete(gradNodes, gradEdges);
    model.sense = repmat('=', size(Aeq, 1), 1);

    params = struct();
    params.Method = 1; % Dual simplex    
    params.OptimalityTol = 1e-4;    
%     params.OutputFlag = 0;
    result = gurobi(model, params);
            
    [xNode, xEdge] = unpackOvercomplete(result.x, L, N, nEdges);
        
    ee = result.objval;        

end

