% script for generating threshold using shuffled data

[~,savename,~] = fileparts(pwd);
savepath = ['/vega/brain/users/sh3276/results/luis/models/shuffled/' savename '.mat'];

% merge_all_models;
load('./results/model_collection.mat');

% initialize
num_shuffle = length(model_collection.models);
theta = model_collection.models{1}.theta;
num_node = size(theta.edge_potentials,1);
node_pot_shuffled = zeros(size(theta.node_potentials,1),num_shuffle);
edge_pot_shuffled = zeros(size(theta.edge_potentials,1),...
    size(theta.edge_potentials,2),num_shuffle);

% collect shuffled data - convert from log scale to real scale
for i = 1:num_shuffle
    grph = model_collection.models{i}.structure;
    theta = model_collection.models{i}.theta;
    node_pot_shuffled(:,i) = theta.node_potentials;
    ep = theta.edge_potentials;
    ep(grph==0) = NaN;
    edge_pot_shuffled(:,:,i) = ep;
    
end

save(savepath,'node_pot_shuffled','edge_pot_shuffled','-v7.3');