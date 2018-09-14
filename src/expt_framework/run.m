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
        fprintf('Loading existing %s', fsave);
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
        x_train_base = X(1:floor(params.split*sample_count),:);
        x_test_base = X((floor(params.split*sample_count)+1):sample_count,:);


        %% Prep data

        x_train = add_lookback_nodes(x_train_base, params.time_span);
        x_test = add_lookback_nodes(x_test_base, params.time_span);
        stim_count = size(params.stimuli, 2);
        % Append any stimulus nodes
        if stim_count > 0
            assert(sample_count == size(params.stimuli, 1), ...
                   'Stimuli and neuron data must have same number of samples.')
            x_train = [x_train params.stimuli(1:floor(params.split*sample_count),:)];
            x_test = [x_test params.stimuli((floor(params.split*sample_count)+1):sample_count,:)];
        end


        % define allowed edges via variable groups
        % One variable group entry per node. Each group entry is a (possibly
        % empty) list of other node indexes.
        if params.time_span > 1
            num_groups = size(x_train,2);
            variable_groups = cell(1, num_groups);
            base_node_count = size(params.data, 2);
            % Indexes of current timestep nodes
            origidx = 1:base_node_count;
            % Indexes of prev timestep nodes
            dupidx = base_node_count+1:(num_groups - stim_count);

            % Always consider fully connected graph at current timestep
            variable_groups(origidx) = all_but_me(1, base_node_count);

            if strcmp(params.edges, 'full')
                % Fully connect all neuron nodes.
                variable_groups = all_but_me(1, params.time_span * base_node_count);
%                     variable_groups = uint16(1:size(x_train, 2));
            else
                if strcmp(params.edges, 'simple')==0
                    fprintf('Invalid edges parameter. Resorting to default setting "simple".\n');
                    params.edges  = 'simple';
                end
                fprintf('Prohibiting edges between offset nodes.\n');
                % Add edge from every current timestep node to every
                % added previous timestep node.

                % Add half edges from orig nodes to every dup node.
                variable_groups(origidx) = cellfun(@(x) [x dupidx], ...
                    variable_groups(origidx),'UniformOutput',false);

                % Set half edge from every dup edge to every orig node.
                variable_groups(dupidx) = {origidx};
            end

            % Always connect stimulus nodes to all non-stimulus nodes
            % NOTE: Does NOT connect stimulus nodes to each other
            stimidx = (num_groups - stim_count + 1):num_groups;
            % Half edges from neuron nodes to stimulus nodes
            for ii = [origidx dupidx]
                variable_groups{ii} = [variable_groups{ii} stimidx];
            end
            % Half edges from stimulus nodes to neuron nodes
            for ii = stimidx
                variable_groups{ii} = [origidx dupidx];
            end

            if params.no_same_neuron_edges
                fprintf('Removing all edges between same-neuron nodes.\n');
                for ii = origidx
                    this_node_idxs = ii:base_node_count:params.time_span*base_node_count;
                    for jj = this_node_idxs
                        variable_groups{jj} = setdiff(variable_groups{jj}, this_node_idxs);
                    end
                end
            end
        else
            % Fully connect ALL nodes (including stimulus)
%                 variable_groups = uint16(1:size(x_train, 2));
            variable_groups = all_but_me(1, size(x_train, 2));
        end


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
                                                    % 'time_span', params.time_span);
            model_collection.variable_groups = variable_groups;
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
        model.time_span = params.time_span;
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

function [indices] = all_but_me(low, high)
%ALL_BUT_ME Cell array of all sets that omit one integer from range [low, high].

    N = high - low + 1;
    tmp = repmat((low:high)', 1, N);
    tmp = tmp(~eye(size(tmp)));
    tmp = reshape(tmp,N - 1, N)';
    indices = num2cell(tmp, 2)';
end
