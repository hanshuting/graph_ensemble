function [ indexes ] = get_ordered_neighborhood_indexes( adjacency_matrix, node_number )
    indexes = find(adjacency_matrix(node_number,:) ~= 0);
end

