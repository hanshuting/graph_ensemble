function [graph,X] = sample_synthetic_data(params)

    % sets the seed
    rand_seed = 1;
    rng('default'); rng(rand_seed);

    % generates a graph
    [node_potentials, edge_potentials] = rand_ising(params.synthetic_graph_size, ...
        params.synthetic_graph_density,'minimum_edge_weight',params.minimum_edge_weight);
    graph.edge_potentials = params.edgeMultiplier*edge_potentials;
    graph.node_potentials = params.edgeMultiplier*node_potentials;
    edge_potentials = full(edge_potentials) * params.edgeMultiplier;
    node_potentials = full(node_potentials) * params.edgeMultiplier;

    % samples from the graph
    X = gibbs_fast_sampler(node_potentials, edge_potentials, params.synthetic_graph_size,...
        params.sample_count, params.burn_in, params.sample_interval, rand_seed);
    X = logical(X);


end

