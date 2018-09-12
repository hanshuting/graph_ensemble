function [ graph ] = chowliu_plot_graph_comparison(x_train, true_adjacency)
    model = treegmFit(x_train);
    
    ref_graph = plot_graph(true_adjacency ~= 0);
    plot_graph(model.adjmat, 'reference_graph', ref_graph);
    
    plot_graph_comparison(true_adjacency ~= 0, model.adjmat);
end

