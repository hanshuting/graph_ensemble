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
    theta = model_collection.models{i}.theta;
    node_pot_shuffled(:,i) = theta.node_potentials;
    edge_pot_shuffled(:,:,i) = theta.edge_potentials;
end

% find significance level
node_thresh = zeros(size(node_pot_shuffled,1),2);
edge_thresh = zeros(size(edge_pot_shuffled,1),size(edge_pot_shuffled,2),2);
for i = 1:num_node
    ccdist = fitdist(exp(node_pot_shuffled(i,:))','normal');
    node_thresh(i,1) = icdf(ccdist,0.05);
    node_thresh(i,2) = icdf(ccdist,0.95);
    for j = 1:num_node
        ccdist = fitdist(reshape(exp(edge_pot_shuffled(i,j,:)),[],1),'normal');
        edge_thresh(i,j,2) = icdf(ccdist,0.05);
        edge_thresh(i,j,2) = icdf(ccdist,0.95);
    end
end

% convert to log scale
node_thresh = log(node_thresh);
edge_thresh = log(edge_thresh);

save(savepath,'edge_thresh','node_thresh','node_pot_shuffled','edge_pot_shuffled','-v7.3');