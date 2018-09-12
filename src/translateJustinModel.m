function mine = translateJustinModel(justin)

    W = sparse(justin.pairs(:,1), ...
               justin.pairs(:,2), ...
               ones(size(justin.pairs, 1), 1), ...
               justin.nnodes, ...
               justin.nnodes);

    W = W + W';
    
    [Aeq, beq] = makeOvercompleteLocalPolytope(W, justin.nvals);
    
    mine = var2struct(Aeq, beq);

end

