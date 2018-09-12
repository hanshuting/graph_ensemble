function [node_potentials, edge_potentials, logZ_potentials] = get_node_and_edge_potentials(F, G, logZ, edge_list)
% Input: 
%   F: 2xNodeCount matrix, but I don't understand meaning (see paper) 
%   G: 4xEdgeCount matrix, but I don't understand meaning (see paper) 
%   edge_list: matrix with one row per edge and two columns indicating the nodes of 
%       the edges, first column is always smaller than the second 
% Output:
%   node_potentials: column array of node potentials
%   edge_potentials: symmetric matrix with edge potentials

    % Node potentials 
    node_count = size(F,2);
    node_potentials = zeros(node_count,1);
    for i = 1:size(G,2)
        node_1 = edge_list(i,1);
        node_2 = edge_list(i,2);
        node_potentials(node_1) = node_potentials(node_1) + G(2,i) - G(1,i);
        node_potentials(node_2) = node_potentials(node_2) + G(3,i) - G(1,i);
    end
    node_potentials = node_potentials + (F(2,:) - F(1,:))';
    
    % Edge potentials 4+1-2-3 (11+00-01-10)
    % Create a row array where each column has the edge potential of one
    % edge
    edge_potentials_array = G(4,:) + G(1,:) - G(2,:) - G(3,:);
    
    edge_potentials = zeros(node_count);
    edge_count = numel(edge_potentials_array);
    
    % transfer edge potentials to upper triangle
    for i = 1:edge_count
        node_1 = edge_list(i,1);
        node_2 = edge_list(i,2);
        edge_potentials(node_1,node_2) = edge_potentials_array(i);
    end
    
    % make matrix symmetric
    edge_potentials = edge_potentials + edge_potentials';
    
    % Update logZ
    % logZ = logZ' - summation of constant terms while converting to node
    % and edge potentials
    logZ_potentials = logZ - sum(F(1,:)) - sum(G(1,:));
end