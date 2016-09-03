function [gcc] = global_cluster_coeff(graph)
% this function calculates the global clustering coefficinet given graph

num_node = size(graph,1);
num_tric = 0;
num_ftric = 0;

% make sure the diagnal is zero
graph = graph-diag(diag(graph));

% go over possible 3-cliques
for i = 1:num_node
    cc_node = find(graph(i,:));
    for j = 1:length(cc_node)-1
        for k = j+1:length(cc_node)
            num_tric = num_tric+1;
            if graph(cc_node(j),cc_node(k))==1
                num_ftric = num_ftric+1;
%                 graph(cc_node(j),cc_node(k)) = 0;
            end
        end
    end
end

gcc = num_ftric/num_tric;

end