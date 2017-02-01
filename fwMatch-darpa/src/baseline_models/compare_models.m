function [ ] = compare_models(graph_size, graph_density, edge_potential_multiplier, sample_count, train_split)
    % sets the seed
    rand_seed = 1;
    rng('default'); rng(rand_seed);
    
    % generates a graph
    [node_potentials, edge_potentials] = rand_ising(graph_size, graph_density);
    edge_potentials = full(edge_potentials) * edge_potential_multiplier;

    % samples from the graph
    X = gibbs_fast_sampler(node_potentials, edge_potentials, graph_size, sample_count, 10000, 100, rand_seed);
    X = logical(X);
    x_train = X(1:floor(train_split*sample_count),:);
    x_test = X(floor(train_split*sample_count)+1:sample_count,:);

    fprintf('\nChow-liu test loglike: %f\n', chowliu_loglike(x_train,x_test));
    fprintf('Bernoulli test loglike: %f\n', bernoulli_loglike(x_train,x_test));
end

