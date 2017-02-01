%% Experiment: Parameter Estimation for Real data
%  Description:
%     - Learn multiple structures with various regularizers (lambda) and
%     tolerance values 
%     - Train the parameters with respect to each of the structure learned
%     in the previous step with various values of parameter estimation
%     regularizer
%     - Predict (test) the learned models in the previous step on validation data
%     - Finally, choose the optimal model based on validation likelihood
%
%  Inputs (with example):
%    params.lambda = 1e-10;
%    params.use_built_in_lasso = false;
%
%  Quick run in MATLAB: 
%    In root folder, after loading startup.m, run:
%    runExpt('structure_learning_synthetic', 1)
%  Quick run on command line:
%    \run.sh Parameter_estimation 1

function run(params)

     %Kate: if result is already computed, we are done
     fsave = sprintf('%s/results/%s', params.exptDir, ...
                     params.saveFileName);
     if(exist(fsave))
       %loads 'model_collection'!!!!
       load(fsave);
       disp(model_collection);
       for i = 1:numel(model_collection.models)
         fprintf('model %d:\n', i);
         model = model_collection.models{i};
         if(isfield(model, 'test_likelihood'))
           continue;
         end
         fprintf('Computing training likelihood\n');
         model.train_likelihood = compute_avg_log_likelihood( ...
             model.theta.node_potentials, ...
             model.theta.edge_potentials, ...
             model.theta.logZ, ...
             model_collection.x_train);
         
         %Kate: reverting changes by rvtonge from 11/11/15:
         % Compute test likelihood
         fprintf('Computing test likelihood\n');
         model.test_likelihood = compute_avg_log_likelihood( ...
             model.theta.node_potentials, ...
             model.theta.edge_potentials, ...
             model.theta.logZ, ...
             model_collection.x_test);
         model_collection.models{i} = model;
       end
     else
       
       %% Get the data
       X = params.data;
       sample_count = size(X,1);
       x_train = X(1:floor(params.split*sample_count),:);
       x_test = X((floor(params.split*sample_count)+1):sample_count,:);
       
       %% Instantiate object that runs the algorithm
       if strcmp(params.structure_type, 'loopy')
         model_collection = LoopyModelCollection(x_train, x_test, ...
                                                 's_lambda_count',params.s_lambda_count, ...
                                                 's_lambda_min',params.s_lambda_min, ...
                                                 's_lambda_max',params.s_lambda_max, ...
                                                 'density_count',params.density_count, ...
                                                 'density_min',params.density_min, ...
                                                 'density_max',params.density_max, ...
                                                 'p_lambda_count',params.p_lambda_count, ...
                                                 'p_lambda_min',params.p_lambda_min, ...
                                                 'p_lambda_max',params.p_lambda_max);
         %model_collection.variable_groups = [ones(1, ];
       else
         % when training a tree there is no need for density and structure
         % lambda
         model_collection = LoopyModelCollection(x_train, x_test, ...
                                                 'p_lambda_count',params.p_lambda_count, ...
                                                 'p_lambda_min',params.p_lambda_min, ...
                                                 'p_lambda_max',params.p_lambda_max);
       end
       
       %% Structure Learninig
       fprintf('Structure Learning...\n');
       if strcmp(params.structure_type, 'loopy')
         % learn loopy models
         model_collection = model_collection.do_loopy_structure_learning();
       else
         % learn tree model (learns single structure)
         model_collection = model_collection.do_chowliu_structure_learning();
       end
     
       %% Parameter Estimation
       fprintf('Parameter Estimation...\n');
       disp(model_collection);
       % Training
       model_collection = model_collection.do_parameter_estimation( ...
           'BCFW_max_iterations', params.BCFW_max_iterations, ...
           'compute_true_logZ', params.compute_true_logZ, ...
           'reweight_denominator', params.reweight_denominator);
     end
    
    % SH
    for i = 1:numel(model_collection.models)
         fprintf('model %d:\n', i);
         model = model_collection.models{i};
         if(isfield(model, 'test_likelihood'))
           continue;
         end
         fprintf('Computing training likelihood\n');
         model.train_likelihood = compute_avg_log_likelihood( ...
             model.theta.node_potentials, ...
             model.theta.edge_potentials, ...
             model.theta.logZ, ...
             model_collection.x_train);

         %Kate: reverting changes by rvtonge from 11/11/15:
         % Compute test likelihood
         fprintf('Computing test likelihood\n');
         model.test_likelihood = compute_avg_log_likelihood( ...
             model.theta.node_potentials, ...
             model.theta.edge_potentials, ...
             model.theta.logZ, ...
             model_collection.x_test);
         model_collection.models{i} = model;
    end

    % Saves results
    disp(model_collection);
    save(sprintf('%s/results/%s', params.exptDir, params.saveFileName), 'model_collection');
end

