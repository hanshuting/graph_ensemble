% save shuffled model results
[~,savename,~] = fileparts(pwd);
savepath = '/vega/brain/users/sh3276/results/luis/models/shuffled/';

load('./results/model_collection.mat');

num_models = length(model_collection.models);
F = cell(num_models,1);
G = cell(num_models,1);
node_pot = cell(num_models,1);
edge_pot = cell(num_models,1);
graphs = cell(num_models,1);
logZ = zeros(num_models,1);

for ii = 1:num_models
    F{ii} = model_collection.models{ii}.theta.F;
    G{ii} = model_collection.models{ii}.theta.G;
    node_pot{ii} = model_collection.models{ii}.theta.node_potentials;
    edge_pot{ii} = model_collection.models{ii}.theta.edge_potentials;
    logZ(ii) = model_collection.models{ii}.theta.logZ;
    graphs{ii} = model_collection.models{ii}.structure;
end

save([savepath savename '_fulldata.mat'],'F','G','node_pot','edge_pot',...
    'graphs','logZ','-v7.3');