
[~,savename,~] = fileparts(pwd);
savepath = './results/';

load('./results/model_collection.mat');
[best_model, best_model_index] = getBestModel(model_collection);
graph = best_model.structure;
node_pot = best_model.theta.node_potentials;
edge_pot = best_model.theta.edge_potentials;
G = best_model.theta.G;
F = best_model.theta.F;
logZ = best_model.theta.logZ;
s_lambda = best_model.s_lambda;
p_lambda = best_model.p_lambda;
density = best_model.density;
time_span = model_collection.models{best_model_index}.time_span;

save([savepath savename '_best_model_full.mat'], 'graph', 'node_pot', 'edge_pot', ...
    'G', 'F', 'logZ', 's_lambda', 'p_lambda', 'density', 'time_span', '-v7.3');
