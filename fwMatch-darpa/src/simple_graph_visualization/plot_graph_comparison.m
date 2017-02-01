function plot_graph_comparison( true_adjacency, adjacency, varargin )
%PLOT_GRAPH_COMPARISON plots the graph highlighting extra and missing edges
%in respect to the true graph
%
% Input
%   true_adjacency: adjacency matrix of the true graph
%   adjacency: adjacency matrix of the estimated graph
%   graph (optional): reference graph to keep node positions

    % parse parameters
    parser = inputParser;
    parser.addRequired('true_adjacency', @ismatrix);
    parser.addRequired('adjacency', @ismatrix);
    parser.addOptional('reference_graph', false);
    parser.parse(true_adjacency, adjacency, varargin{:});
    reference_graph = parser.Results.reference_graph;

    true_adjacency = logical(true_adjacency ~= 0);
    adjacency = logical(adjacency ~= 0);
    node_count = length(true_adjacency);
    node_labels = regexp(num2str(linspace(1,node_count,node_count)),'\s+','split');

    % create graph with all edges (from true and predicted graphs)
    graph = biograph(triu(true_adjacency | adjacency),node_labels,'ShowArrows','off','ShowWeights','off','ShowTextInNodes','label');
    set(graph.Edges,'LineWidth',1);
    set(graph.Nodes, 'Shape', 'circle');
    
    % color true positives green
    true_positives = (adjacency & true_adjacency);
    [rows,cols] = find(triu(true_positives));
    node_ids = reshape(strtrim(cellstr(num2str([rows; cols]))), size(rows,1), 2);
    true_positive_edges = getedgesbynodeid(graph, node_ids);
    set(true_positive_edges,'LineColor',[.4 1 .4]);
    
    % color false positives grey
    false_positives = (adjacency - true_adjacency > 0);
    [rows,cols] = find(triu(false_positives));
    if length(rows) > 0
        node_ids = reshape(strtrim(cellstr(num2str([rows; cols]))), size(rows,1), 2);
        false_positive_edges = getedgesbynodeid(graph, node_ids);
        set(false_positive_edges,'LineColor',[.5 .5 .5]);
    end
    
    % color false negatives red
    false_negatives = (adjacency - true_adjacency < 0);
    [rows,cols] = find(triu(false_negatives));
    if length(rows) > 0
        node_ids = reshape(strtrim(cellstr(num2str([rows; cols]))), size(rows,1), 2);
        false_negative_edges = getedgesbynodeid(graph, node_ids);
        set(false_negative_edges,'LineColor',[1 0.2 0.2]);
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

