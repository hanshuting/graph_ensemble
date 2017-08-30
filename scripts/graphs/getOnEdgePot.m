function [G_on] = getOnEdgePot(graph,G)

num_node = size(graph,1);
trilgraph = tril(graph);
num_edge = sum(sum(logical(trilgraph)));
edge_list = zeros(num_edge,2);
[edge_list(:,2),edge_list(:,1)] = find(trilgraph);
G_on = zeros(num_node,num_node);
for i = 1:size(edge_list,1)
    node_1 = edge_list(i,1);
    node_2 = edge_list(i,2);
    G_on(node_1,node_2) = G(4,i);%-G(3,i)-G(2,i);
end

end