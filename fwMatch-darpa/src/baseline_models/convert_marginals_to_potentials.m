function [node_pot,edge_pot] = convert_marginals_to_potentials(x_train, p_lambda, tree_structure, node_marginals, edge_marginals, edge_list)

overcomplete_struct = samples_to_overcomplete(x_train, tree_structure);

mu = reshape(node_marginals, [size(x_train,2), 1]);
TN = [1-mu mu];

TE = [];
for i=1:length(edge_marginals)
    if edge_list(i,1) > edge_list(i,2)
        TE = [TE;reshape(edge_marginals(:,:,i)',[1 4])];
        edge_list(i,:) = [edge_list(i,2) edge_list(i,1)];
    else
        TE = [TE;reshape(edge_marginals(:,:,i)',[1 4])];
    end       
end
[~, edge_indices] = sortrows(edge_list);
TE = TE(edge_indices,:);

F = -(1/p_lambda) * (overcomplete_struct.Ut*(repmat(TN, size(x_train,1),1) - overcomplete_struct.YN)).';
G = -(1/p_lambda) * (overcomplete_struct.Vt*(repmat(TE, size(x_train,1),1) - overcomplete_struct.YE)).';

[node_pot, edge_pot, ~] = get_node_and_edge_potentials...
            (F, G, 0, overcomplete_struct.edges{1}');
end