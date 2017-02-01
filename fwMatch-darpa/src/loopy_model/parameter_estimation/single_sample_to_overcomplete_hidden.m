function [y_node, y_edge] = single_sample_to_overcomplete_hidden(sample, label_count, edges)
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
        y_node(:, label) = (sample.^(label-1)) .* ((1-sample).^(2-label));
    end
            
    % create overcomplete representation for edges
    y_edge = zeros(edge_count, label_count^2);
    [x_i,x_j] = ind2sub([label_count, label_count],1:(label_count^2));
    for edge_index = 1:edge_count
        i = edges(edge_index,1);
        j = edges(edge_index,2);
        
        y_i = sample(i);
        y_j = sample(j);
        
        y_edge(edge_index, :) = (y_i.^(x_i-1)) .* ((1-y_i).^(2-x_i)) .*...
            (y_j.^(x_j-1)) .* ((1-y_j).^(2-x_j));        
    end
    
    % Check our work
    assert(all(abs(sum(y_node, 2)-1) <= eps(1)), 'y_node not correctly labeled.');
    assert(all(abs(sum(y_edge, 2)-1) <= eps(1)), 'y_edge not correctly labeled.');
end
