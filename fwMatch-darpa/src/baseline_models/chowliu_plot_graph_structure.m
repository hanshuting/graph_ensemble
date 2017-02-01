function [ graph ] = chowliu_plot_graph_structure(x_train, node_names)
    model = treegmFit(x_train);
    
    if nargin > 1
        plot_graph(model.adjmat, 'node_names', node_names);
    else
        plot_graph(model.adjmat);
    end
end

