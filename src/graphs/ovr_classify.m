
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
EXPT = 'm21_sample_d2_d3';
EE = {'vis_01','vis_02'};
LPATH = '/vega/brain/users/sh3276/src/fwMatch-darpa/expt';
SPATH = '/vega/brain/users/sh3276/results/luis';

% test info
TEST_EE = {'po_vis_01','po_vis_02'};
DPATH = '/vega/brain/users/sh3276/data/luis';

%% calculate LL with each model
for nn = 1:length(TEST_EE)
    
    % load test data
    fload = sprintf('%s/%s_%s.mat',DPATH,EXPT,TEST_EE{nn});
    load(fload);
    fprintf('loaded test file: %s\n',fload);
    samples_N = size(data,1);
    if(islogical(data) || length(unique(data(:))) == 2)
        binary_model = 1;
    else
        binary_model = 0;
    end

    LLs = zeros(size(data,1),length(EE));

    for ii = 1:length(EE)
    
        % load model
        fload = sprintf('%s/%s_%s_%s/results/model_collection.mat',LPATH,EXPT,...
            EE{ii},MODEL_TYPE);
        load(fload);
        fprintf('loaded model %s\n',fload);
    
        best_model = getBestModel(model_collection);
        inferred_model = best_model.inference_model;
        clear model_collection best_model;
    
        node_pot = inferred_model.theta.node_potentials;
        edge_pot = inferred_model.theta.edge_potentials;
        if(isfield(inferred_model.theta, 'true_logZ') )
            logZ = inferred_model.theta.true_logZ;
        else
            logZ = inferred_model.theta.logZ;
        end
        clear inferred_model;
    
        % % threshold potentials
        if if_shuffle
            fload = sprintf('%s/models/shuffled/shuffled_%s_%s_%s.mat',SPATH,...
                EXPT,EE{ii},MODEL_TYPE);
            load(fload);
            fprintf('loaded shuffled threshold %s\n',fload);
            
            node_pot = node_pot.*double(node_pot>node_thresh);
            edge_pot = edge_pot.*double(edge_pot>edge_thresh);
        end
        
        for s = 1:samples_N
            if(binary_model == 1)
                LLs(s,ii) = compute_avg_log_likelihood(node_pot,edge_pot,logZ,data(s,:));
            else
                LLs(s,ii) = compute_avg_log_likelihood_hidden(node_pot,edge_pot,logZ,data(s,:));
            end
        %   if(mod(s,500) == 0)
        %     fprintf('done test frame %d/%d\n', s, samples_N);
        end
    
    end
    
    % plot
    figure;set(gcf,'color','w');
    plot(1:samples_N,LLs');
    xlabel('frame');ylabel('LL');
    l = legend(EE{:}); set(l,'Interpreter','none');
    title([TEST_EE{nn} '_inference'],'Interpreter', 'none');
    fsave = sprintf('%s/infer/LLs_%s_%s.fig',SPATH,EXPT,TEST_EE{nn});
    saveas(gcf,fsave);
    
    % save
    fsave = sprintf('%s/infer/LLs_%s_%s.mat',SPATH,EXPT,TEST_EE{nn});
    save(fsave, 'LLs', 'EXPT');
    fprintf('SAVED: %s\n', fsave);

end
