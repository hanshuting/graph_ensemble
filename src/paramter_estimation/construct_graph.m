function [graphUT,graphComplete] = construct_graph(edge_weights,num_nodes)

graph = zeros(num_nodes);

for i=1:num_nodes-1
    for j=i+1:num_nodes
        graph(i,j) = edge_weights(getEdgeId(num_nodes,[i,j]));
    end
end
graphUT = graph;
graphComplete = graph + graph';

end
