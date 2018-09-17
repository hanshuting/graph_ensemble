function [G_on] = getOnEdgePot(graph,G)
%GETONEDGEPOT Return just the (1, 1) edge potentials for all node pairs.
%
%   Returns
%       G_on: NxN matrix, but only upper triangle is populated.
num_node = size(graph,1);
trilgraph = triu(graph); % SH 5/5/2018
num_edge = sum(sum(logical(trilgraph)));
edge_list = zeros(num_edge,2);
[edge_list(:,2),edge_list(:,1)] = find(trilgraph);
G_on = zeros(num_node,num_node);
for i = 1:num_edge
    node_1 = edge_list(i,1);
    node_2 = edge_list(i,2);
    G_on(node_1,node_2) = G(4,i);%-G(3,i)-G(2,i);
end

end