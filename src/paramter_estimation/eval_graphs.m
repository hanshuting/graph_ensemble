function eval = eval_graphs(true_graph,estimated_graph)

eval.TP  = sum(sum(triu(true_graph) & triu(estimated_graph)));
eval.FP  = sum(sum(~triu(true_graph) & triu(estimated_graph)));
eval.FN  = sum(sum(triu(true_graph) & ~triu(estimated_graph)));

end

