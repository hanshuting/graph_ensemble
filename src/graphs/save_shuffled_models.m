% save shuffled model results
[~,savename,~] = fileparts(pwd);
savepath = './results/';

load('./results/model_collection.mat');

num_models = length(model_collection.models);
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
    F{ii} = model_collection.models{ii}.theta.F;
    G{ii} = model_collection.models{ii}.theta.G;
    node_pot{ii} = model_collection.models{ii}.theta.node_potentials;
    edge_pot{ii} = model_collection.models{ii}.theta.edge_potentials;
    logZ(ii) = model_collection.models{ii}.theta.logZ;
    graphs{ii} = model_collection.models{ii}.structure;
    s_lambda(ii) = model_collection.models{ii}.s_lambda;
    p_lambda(ii) = model_collection.models{ii}.p_lambda;
    density(ii) = model_collection.models{ii}.density;
    time_span(ii) = model_collection.models{ii}.time_span;
end

save([savepath savename '_fulldata.mat'],'F','G','node_pot','edge_pot',...
    'graphs','logZ', 's_lambda', 'p_lambda', 'density','time_span','-v7.3');
