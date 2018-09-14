function [subgraph] = high_score_subgraph(graph,G)

num_iter = 100;

% get edge list
graph = tril(graph);
num_node = size(graph,1);
num_edge = sum(sum(logical(graph)));
edge_list = zeros(num_edge,2);
[edge_list(:,2),edge_list(:,1)] = find(graph);

% get synchronous term
G_on = zeros(num_node,num_node);
for i = 1:size(G,2)
    node_1 = edge_list(i,1);
    node_2 = edge_list(i,2);
    G_on(node_1,node_2) = G(4,i);
%     G_on(node_1,node_2) = G(4,i)-G(3,i)-G(2,i);
end

% loss function
edge_h = @(p,edge_p) -sum(sum(edge_p(logical(p),logical(p))));
loss = @(x) edge_h(x,G_on);

% simulated annealing
opt.Verbosity = 0;
sol = cell(num_iter,1);
fval = zeros(num_iter,1);
for i = 1:num_iter
    init_guess = rand(1,num_node)>0.7;
    [minimum,fval(i)] = anneal(loss,init_guess,opt);
    sol{i} = find(minimum);
end

% find the best solution
min_indx = find(fval==min(fval));
min_sol = sol(min_indx);
if length(min_indx)>1
%     sz_sol = cellfun('length',min_sol);
%     best_indx = find(sz_sol==min(sz_sol));
%     if length(best_indx)>1
        core_indx = min_sol{1};
        for i = 1:length(min_indx)
            core_indx = intersect(core_indx,min_sol{i});
%             core_indx = unique([core_indx,min_sol{i}]);
        end
        subgraph = core_indx;
%     else
%         subgraph = min_sol{best_indx};
%     end
else
    subgraph = min_sol{1};
end

end