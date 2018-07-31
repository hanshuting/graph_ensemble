function save_and_merge_shuffled_models(results_folder)
%SAVE_AND_MERGE_SHUFFLED_MODELS
merged_models = merge_all_results(results_folder);

num_models = length(merged_models.models);
F = cell(num_models,1);
G = cell(num_models,1);
node_pot = cell(num_models,1);
edge_pot = cell(num_models,1);
graphs = cell(num_models,1);
logZ = zeros(num_models,1);
s_lambda = zeros(num_models,1);
p_lambda = zeros(num_models,1);
density = zeros(num_models,1);
time_span = zeros(num_models,1);

for ii = 1:num_models
    F{ii} = merged_models.models{ii}.theta.F;
    G{ii} = merged_models.models{ii}.theta.G;
    node_pot{ii} = merged_models.models{ii}.theta.node_potentials;
    edge_pot{ii} = merged_models.models{ii}.theta.edge_potentials;
    logZ(ii) = merged_models.models{ii}.theta.logZ;
    graphs{ii} = merged_models.models{ii}.structure;
    s_lambda(ii) = merged_models.models{ii}.s_lambda;
    p_lambda(ii) = merged_models.models{ii}.p_lambda;
    density(ii) = merged_models.models{ii}.density;
    time_span(ii) = merged_models.models{ii}.time_span;
end

save(fullfile(results_folder, 'fulldata.mat'),'F','G','node_pot','edge_pot',...
    'graphs','logZ', 's_lambda', 'p_lambda', 'density','time_span','-v7.3');

end
