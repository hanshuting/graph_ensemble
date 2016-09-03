% extract ensemble data for training

rng(1000);
expt_name = {'m21_d2_vis'};
ee = {{'01_high'}};
model_ee = {'loopy_mem_core'};
%{'svd_core','cc_mem_core','loopy_mem_core','cc_cent_core',...
%    'loopy_cent_core'};
%    ,'svd_noncore','cc_mem_noncore','loopy_mem_noncore',...
%    'cc_cent_noncore','loopy_cent_noncore'};
test_ee = {{'01_high','02_high','off_high','spon_high'}};
k = 3;
num_shuff = 100;

%% calculate LL
for n = 1:length(expt_name)
    
    expt_ee = ee{n};
    expt_test = test_ee{n};
    
    for e = 1:length(expt_ee)
        for m = 1:length(model_ee)
            
            fprintf('running experiment %s_%s_%s_k_%u...\n',expt_name{n},...
                expt_ee{e},model_ee{m},k);
            
            % set path
            data_path = ['/vega/brain/users/sh3276/data/luis/']; 
            model_path = ['/vega/brain/users/sh3276/src/fwMatch-darpa/expt/'];
            shuff_path = ['/vega/brain/users/sh3276/src/fwMatch-darpa/expt/'];
            save_path = ['/vega/brain/users/sh3276/results/luis/infer/' ...
                expt_name{n} '/'];
            
            load([model_path expt_name{n} '_' expt_ee{e} '_' model_ee{m} ...
                '_train_k_' num2str(k) '_loopy/results/model_collection.mat']);
            best_model = getBestModel(model_collection);
            inferred_model = best_model.inference_model;
            node_pot = inferred_model.theta.node_potentials;
            edge_pot = inferred_model.theta.edge_potentials;
            if(isfield(inferred_model.theta, 'true_logZ') )
                logZ = inferred_model.theta.true_logZ;
            else
                logZ = inferred_model.theta.logZ;
            end
            clear inferred_model model_collection;
            
            % LL with CRF model and with shuffled model
            LL = cell(length(expt_test),1);
            LL_shuff = cell(length(expt_test),1);
            for t = 1:length(expt_test)
                
                %---------- CRF ----------%
                if ~strcmp(expt_ee{e},expt_test{t})
                    load([data_path expt_name{n} '_' expt_ee{e} '_' ...
                        model_ee{m} '_' expt_test{t} '_k_' num2str(k) '.mat']);
                else
                    load([data_path expt_name{n} '_' expt_ee{e} '_' ...
                        model_ee{m} '_test_k_' num2str(k) '.mat']);
                end
                
                num_sample = size(data,1);
                test_LL = zeros(num_sample,1);
                
                % works only for binary model now
                for s = 1:num_sample
                    test_LL(s) = compute_avg_log_likelihood(node_pot,...
                        edge_pot,logZ,data(s,:));
                end
                
                LL{t} = test_LL;
                
                %---------- shuffled ----------%
                load([shuff_path 'shuffled_' expt_name{n} '_' expt_ee{e}...
                    '_' model_ee{m} '_train_k_' num2str(k) ...
                    '_loopy/results/model_collection.mat']);
                test_LL_shuff = zeros(num_sample,num_shuff);
                for shuff = 1:num_shuff
                    single_model =  SingleLoopyModel( ...
                        model_collection.x_train, model_collection.x_test, ...
                        model_collection.models{shuff}, ...
                        model_collection.variable_names);
                    node_pot = single_model.theta.node_potentials;
                    edge_pot = single_model.theta.edge_potentials;
                    logZ = single_model.theta.logZ;
                    for s = 1:num_sample
                        test_LL_shuff(s,shuff) = compute_avg_log_likelihood...
                            (node_pot,edge_pot,logZ,data(s,:));
                    end
                end
                LL_shuff{t} = test_LL_shuff;
                
            end
            
            save([save_path expt_name{n} '_' expt_ee{e} '_' model_ee{m} ...
                '_k_' num2str(k) '_LL.mat'],'LL','LL_shuff','expt_test');
            
        end
    end
end
