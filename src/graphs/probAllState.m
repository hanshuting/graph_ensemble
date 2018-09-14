rng(1000);
MODEL_TYPE = 'loopy';
EXPT = 'm21_d2_vis';
EE = {'01_high'};
LPATH = '/vega/brain/users/sh3276/src/fwMatch-darpa/expt';
SPATH = '/vega/brain/users/sh3276/results/luis/infer';

% test info
num_test = 1e5;
DPATH = '/vega/brain/users/sh3276/data/luis';

for i = 1:length(EE)
    
    % load model
    fload = sprintf('%s/%s_%s_%s/results/model_collection.mat',LPATH,EXPT,...
        EE{i},MODEL_TYPE);
    load(fload);
    fprintf('loaded model %s\n',fload);

    best_model = getBestModel(model_collection);
    inferred_model = best_model.inference_model;
    clear model_collection best_model;

    graph = inferred_model.structure;
    num_node = size(graph,1);
    node_pot = inferred_model.theta.node_potentials;
    edge_pot = inferred_model.theta.edge_potentials;
    F = inferred_model.theta.F;
    G = inferred_model.theta.G;
    if(isfield(inferred_model.theta, 'true_logZ') )
        logZ = inferred_model.theta.true_logZ;
    else
        logZ = inferred_model.theta.logZ;
    end
    clear inferred_model;
    
%     LL = zeros(num_test,1);
%     LL_state = zeros(num_node,num_node,num_test,'uint8');
    LL = [];
    LL_state = [];
    num_active = 0:20;
    for nn = 1:length(num_active)
        state_vec = randi([0 1],num_node,num_node);
        LL(j) = compute_sample_likelihood_overcomplete(F,G,logZ,state_vec,graph);
        LL_state(:,:,j) = state_vec;
    end
    
    LL_state = logical(LL_state);
    save(sprintf('%s/%s_%s_LL_states.mat',SPATH,EXPT,EE{i}),'LL','LL_state','-v7.3');
    
end
