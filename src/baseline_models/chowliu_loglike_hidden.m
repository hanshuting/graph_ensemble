function best_test_loglike = chowliu_loglike_hidden(x_train, x_test, p_lambda_sequence)
    model = treegmFitBinary(x_train, [0 1]);
    
    [CPT,edge_list] = ConvertChowLiuModelToPsi(model);
    [~, node_marginals, ~, edge_marginals] = solveConditionalDAITree(CPT,edge_list,model.Nnodes);
    
    best_test_loglike = -Inf;
    for p_lambda = p_lambda_sequence
        [node_pot,edge_pot] = convert_marginals_to_potentials...
            (x_train, p_lambda, model.adjmat, node_marginals, edge_marginals, edge_list);
        [~,~,~,~,~,~,true_logZ] = run_junction_tree(node_pot, edge_pot, 'verbose', true);
    
        best_test_loglike = max(compute_avg_log_likelihood_hidden(node_pot, edge_pot, true_logZ, x_test)...
            , best_test_loglike);
    end
end

