
% add required path
if(isempty(which('LoopyModelCollection')))
  addpath(genpath('/vega/brain/users/sh3276/src/fwMatch/'));
  exptPath = pwd;
  cd('/vega/brain/users/sh3276/src/fwMatch/');
  startup;
  cd(exptPath);  
end
dbclear all

%% set parameters
if_shuffle = 0;

% model info
MODEL_TYPE = 'loopy';
EXPT = 'm21_d2_vis';
EE = {'01_high_tf','02_high_tf'};
LPATH = '/vega/brain/users/sh3276/src/fwMatch-darpa/expt';
SPATH = '/vega/brain/users/sh3276/results/luis/infer';

% test info
TEST_EE = {'all_high_tf'};
DPATH = '/vega/brain/users/sh3276/data/luis';

%% calculate LL with each model
for nn = 1:length(TEST_EE)
    
    % load test data
    fload = sprintf('%s/%s_%s.mat',DPATH,EXPT,TEST_EE{nn});
    load(fload);
    fprintf('loaded test file: %s\n',fload);
    samples_N = size(data,1);
    num_spikes = sum(data,2);
    if(islogical(data) || length(unique(data(:))) == 2)
        binary_model = 1;
    else
        binary_model = 0;
    end
    
    for ii = 1:length(EE)
    
        % load model
        fload = sprintf('%s/%s_%s_%s/results/model_collection.mat',LPATH,EXPT,...
            EE{ii},MODEL_TYPE);
        load(fload);
        fprintf('loaded model %s\n',fload);
        
        num_model = length(model_collection.models);
        LLs = zeros(num_model,samples_N);
        params = zeros(num_model,3);
        for jj = 1:num_model
            
            single_model = SingleLoopyModel( ...
                model_collection.x_train, model_collection.x_test, ...
                model_collection.models{jj}, ...
                model_collection.variable_names);
            inferred_model = single_model.inference_model;
            node_pot = inferred_model.theta.node_potentials;
            edge_pot = inferred_model.theta.edge_potentials;
            
            if(isfield(inferred_model.theta, 'true_logZ') )
                logZ = inferred_model.theta.true_logZ;
            else
                logZ = inferred_model.theta.logZ;
            end
            
            % parameters
            params(jj,1) = inferred_model.s_lambda;
            params(jj,2) = inferred_model.p_lambda;
            params(jj,3) = inferred_model.density;
            
            % likelihood
            for s = 1:samples_N
                LLs(jj,s) = compute_avg_log_likelihood(node_pot,edge_pot,logZ,data(s,:));
            end
            
            % plot
%             h = figure;set(gcf,'color','w');
%             plot(exp(LLs(jj,:)),'b');hold on;plot(num_spikes/max(num_spikes),'r');
%             title(['model ' num2str(jj)]);
            
        end
        
        save(sprintf('%s/%s_%s_infer_%s_all_model_LL.mat',SPATH,EXPT,EE{ii},...
            TEST_EE{nn}),'LLs','params');
        
    end

end
