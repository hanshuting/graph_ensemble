function [graph] = plot_graph(adjacency_matrix, varargin)
%PLOT_GRAPH Plots a small graph (up to 50 nodes perhaps), given its adjacency matrix
%
%   Examples: 
%       plot_graph(adjacency_matrix)
%       plot_graph(adjacency_matrix, 'node_names', {'A', 'B'})
%       plot_graph(adjacency_matrix, 'color_edges', true)
%       plot_graph(adjacency_matrix, 'show_weights', true)
%       plot_graph(adjacency_matrix, 'reference_graph', graph) % keep node positions
    
    % parse input arguments
    parser = inputParser;
    parser.addRequired('adjacency_matrix', @ismatrix);
    parser.addParamValue('node_names', {}, @iscellstr);
    parser.addParamValue('color_edges', false, @islogical);
    parser.addParamValue('show_weights', false, @islogical);
    parser.addParamValue('reference_graph', false);
    parser.parse(adjacency_matrix, varargin{:});
    node_names = parser.Results.node_names;
    color_edges = parser.Results.color_edges;
    show_weights = parser.Results.show_weights;
    reference_graph = parser.Results.reference_graph;
    
    if length(node_names) == 0
        node_names = regexp(num2str(linspace(1,length(adjacency_matrix),length(adjacency_matrix))),'\s+','split');
    end
    
    if show_weights
        show_weights = 'on';
    else
        show_weights = 'off';
    end
    
    graph = biograph(tril(adjacency_matrix), node_names, ...
                     'ShowArrows', 'off', ...
                     'ShowWeights', show_weights, ...
                     'ShowTextInNodes','label', ...
                     'NodeAutoSize', 'on');
    set(graph.Edges,'LineWidth',1);
    %set(graph.Nodes, 'Shape', 'circle');
    
    if color_edges
        max_edge_potential = max(max(adjacency_matrix));
        min_edge_potential = min(min(adjacency_matrix));
        
        for i = 1:length(graph.edges)
            [node_a_id node_b_id] = strtok(graph.edges(i).ID, ' -> ');
            node_a_id = strtrim(node_a_id);
            node_b_id = node_b_id(5:end);
            node_b_id = strtrim(node_b_id);
            
            [ans node_a_index] = ismember(node_a_id, node_names);
            [ans node_b_index] = ismember(node_b_id, node_names);
            edge_potential = adjacency_matrix(node_a_index, node_b_index);
            
            if edge_potential > 0
                relative_edge_potential = edge_potential / max_edge_potential;
                set(graph.edges(i),'LineColor',[(1-relative_edge_potential)*.95 (1-relative_edge_potential)*.95 1]);
            elseif edge_potential < 0
                relative_edge_potential = edge_potential / min_edge_potential;
                set(graph.edges(i),'LineColor',[1 (1-relative_edge_potential)*.95 (1-relative_edge_potential)*.95]);
            end
        end
    end
    
    for i = 1:length(graph.nodes)
        graph.nodes(i).FontSize = 36;
    end
    
    % sets node positions according to reference graph
    if reference_graph ~= false
        % computes initial low energy node positions
        dolayout(reference_graph);
        dolayout(graph);
        
        % copy node positions from reference_graph to graph
        for i = 1:length(graph.nodes)
           graph.nodes(i).position = reference_graph.nodes(i).position;
        end
        
        % recompute edges only
        dolayout(graph, 'Paths', true);
    end
    
    h = view(graph);
end

