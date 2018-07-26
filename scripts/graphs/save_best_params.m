function save_best_params(results_folder)
%SAVE_BEST_PARAMS Find and save the best model and parameters.
%   Writes a single compiled mat file of the result files to disk. Also
%   writes the best model to disk in a mat file, including its parameters.
merged_models = merge_all_results(results_folder);

[best_model, best_model_index] = getBestModel(merged_models);
graph = best_model.structure;
node_pot = best_model.theta.node_potentials;
edge_pot = best_model.theta.edge_potentials;
G = best_model.theta.G;
F = best_model.theta.F;
logZ = best_model.theta.logZ;
s_lambda = best_model.s_lambda;
p_lambda = best_model.p_lambda;
density = best_model.density;
time_span = merged_models.models{best_model_index}.time_span;

save([results_folder 'best_model_full.mat'], 'graph', 'node_pot', 'edge_pot', ...
    'G', 'F', 'logZ', 's_lambda', 'p_lambda', 'density', 'time_span');

save([results_folder 'best_parameters.txt'], 's_lambda', 'density', 'p_lambda', ...
    'time_span', '-ascii', '-double');

end

