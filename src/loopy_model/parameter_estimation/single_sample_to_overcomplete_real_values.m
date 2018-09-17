function [y_node, y_edge] = single_sample_to_overcomplete_real_values(sample, label_count, edges)
% overcompletesample  Convert node sample to overcomplete indicators
%
%   Input:
%       sample: column vector with labels {0,..,label_count} for each variable/node
%       label_count: number of possible labels
%       edges: matrix with one row per each, each row has two node numbers
%           to represent the edge. The smaller node is always in the first
%           column.
%
%   Output:
%       y_node: matrix with one row per variable. Each row is an indicator
%           vector for the label
%       y_edge: matrix with one row per edge. Each row is an indicator
%           vector corresponding to the combination of labels of the two
%           nodes connected through that edge

    assert(isvector(sample), 'sample must be vector; did you pass in matrix?');

    node_count = length(sample);
    edge_count = size(edges, 1);
    
    % create overcomplete representation for nodes
    y_node = zeros(node_count, label_count);
    for label = 1:label_count
        y_node(sample == label-1, label) = 1;
    end
    % Accomodate real value
    y_node(sample == 0.5, :) = 0.5;
            
    % create overcomplete representation for edges
    y_edge = zeros(edge_count, label_count^2);
    mapping_2D_to_1D = reshape(1:(label_count^2), label_count, label_count);
    for edge_index = 1:edge_count
        node_1 = edges(edge_index,1);
        node_2 = edges(edge_index,2);
        
        node_1_label = sample(node_1);
        node_2_label = sample(node_2);
        
        if isequal(fix(node_1_label),node_1_label) && isequal(fix(node_2_label),node_2_label) 
            index = mapping_2D_to_1D(node_1_label+1, node_2_label+1);
            y_edge(edge_index, index) = 1;
        elseif ~isequal(fix(node_1_label),node_1_label) && ~isequal(fix(node_2_label),node_2_label) 
            y_edge(edge_index, :) = 0.25;
        elseif ~isequal(fix(node_1_label),node_1_label)
            index1 = mapping_2D_to_1D(1, node_2_label+1);
            index2 = mapping_2D_to_1D(2, node_2_label+1);
            y_edge(edge_index, index1) = 1;
            y_edge(edge_index, index2) = 1;
            y_edge(edge_index, :) = y_edge(edge_index, :) / 2;
        else
            index1 = mapping_2D_to_1D(node_1_label+1, 1);
            index2 = mapping_2D_to_1D(node_1_label+1, 2);
            y_edge(edge_index, index1) = 1;
            y_edge(edge_index, index2) = 1;
            y_edge(edge_index, :) = y_edge(edge_index, :) / 2;
        end    
    end
    
    % Check our work
    assert(all(sum(y_node, 2) == 1), 'y_node not correctly labeled.');
    assert(all(sum(y_edge, 2) == 1), 'y_edge not correctly labeled.');
end
