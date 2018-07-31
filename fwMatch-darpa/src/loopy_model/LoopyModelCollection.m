classdef LoopyModelCollection

    properties
        % training and test samples
        x_train, x_test;

        % cell array with strings for each node
        variable_names;

        % TODO
        variable_groups;

        % configuration for training
        s_lambda_sequence, p_lambda_sequence, density_sequence;

        % Struct array of estimated parameters of each model
        % Fields:
        %   theta: struct with
        %       F, G (see paper to understand),
        %       node_potentials (column vector),
        %       edge_potentials (symmetric matrix),
        %       true_logZ (computed with JTA)
        %       logZ (computed with Bethe approx.)
        %   s_lambda: regularization of structure learning
        %   p_lambda: regularization of parameter learning
        %   density: target density of the structure
        %   structure: NxN binary adjacency matrix
        %   max_degree: maximum number of connection to a node
        %   train_likelihood: avg. (per sample) likelihood of training set
        %   test_likelihood: avg. (per sample) likelihood of test set
        %   true_test_likelihood: avg. (per sample) likelihood of test set
        %   true_node_marginals: marginals stored when we run JTA to get
        %       true_logZ. These marginals are used to mix the loopy model
        %       with the temporal models
        models;

        % boolean indicator to know if true logZ is to be computed in this
        % run (sometimes we can't because the model is too dense and it
        % would take too long)
        computed_true_logZ;

        % boolean indicator to know whether to use hidden model (for real values)
        % or use binary model
        hidden_model;
    end

    methods
        function self = LoopyModelCollection(x_train, x_test, varargin)
            % Inputs (Required):
            %   x_train: training_samples, one sample per row, logical values
            %   x_test: test_samples, one sample per row, logical values
            % Inputs (Optional: training configuration):
            %   s_lambda_count: number of structure learning regularizers to try
            %   s_lambda_min: minimum lambda value to try
            %   s_lambda_max: maximum lambda value to try
            %   density_count: number of structure learning target densities to try
            %   density_min: minimum density value to try
            %   density_max: maximum density value to try
            %   p_lambda_count: number of parameter learning regularizers to try
            %   p_lambda_min: minimum lambda value to try
            %   p_lambda_max: maximum lambda value to try
            % Inputs (Optional: others)
            %   variable_names: names of the variables (columns of samples)

            % parse input creating an inputParser object
            parser = LoopyModelCollection.parse_and_validate_input(x_train, x_test, varargin{:});

            % initialize training and test sets
            self.x_train = parser.Results.x_train;
            self.x_test = parser.Results.x_test;
            self.variable_names = parser.Results.variable_names;

            % structure lambdas are samples in a logspace
            s_lambda_count = parser.Results.s_lambda_count;
            s_lambda_min_exp = log10(parser.Results.s_lambda_min);
            s_lambda_max_exp = log10(parser.Results.s_lambda_max);
            self.s_lambda_sequence = logspace(s_lambda_min_exp, s_lambda_max_exp, s_lambda_count);

            % parameter lambdas are samples in a logspace
            p_lambda_count = parser.Results.p_lambda_count;
            p_lambda_min_exp = log10(parser.Results.p_lambda_min);
            p_lambda_max_exp = log10(parser.Results.p_lambda_max);
            self.p_lambda_sequence = logspace(p_lambda_min_exp, p_lambda_max_exp, p_lambda_count);

            % density are samples in a linear space
            density_count = parser.Results.density_count;
            density_min = parser.Results.density_min;
            density_max = parser.Results.density_max;
            self.density_sequence = linspace(density_min, density_max, density_count);

            self.models = {};
        end

        % runs the structure learning algorithm with all regularization
        % parameters in s_lambda_sequence and density_sequence
        function self = do_loopy_structure_learning(self)
            fprintf('Structure learning using Lasso Logistic Regression\n');

            % for every s_lambda and density let's learn a structure
            for i = 1:numel(self.s_lambda_sequence)
                % a single call learns structures for all densities
                learned_structures = learn_structures_by_density(self.x_train, ...
                    self.s_lambda_sequence(i), self.density_sequence, ...
                    self.variable_groups);

                % we will initialize structures for all combination of
                % s_lambda, density and p_lambda now. This means we will
                % replicate the same structure |p_lambda_sequence| times
                for j = 1:numel(self.density_sequence)
                    for k = 1:numel(self.p_lambda_sequence)
                        model = struct();
                        model.s_lambda = self.s_lambda_sequence(i);
                        model.density = self.density_sequence(j);
                        model.p_lambda = self.p_lambda_sequence(k);
                        model.structure = learned_structures{j};
                        model.max_degree = max(sum(model.structure));
                        model.median_degree = median(sum(model.structure));
                        model.mean_degree = mean(sum(model.structure));
                        model.rms_degree = rms(sum(model.structure));
                        model.pending_parameter_estimation = true;

                        self.models{end+1} = model;
                    end
                end
            end
        end

        % runs the structure learning algorithm using chowliu. This means
        % there is going to be only one structure
        function self = do_chowliu_structure_learning(self)
            fprintf('Structure learning using Chow-Liu\n');

            % run chowliu once to get adjacecy matrix
            chowliu_output = treegmFitBinary(self.x_train, [0, 1]);
            chowliu_structure = chowliu_output.adjmat;

            for k = 1:numel(self.p_lambda_sequence)
                model = struct();
                model.p_lambda = self.p_lambda_sequence(k);
                model.s_lambda = 0;
                model.density = 0;
                model.structure = chowliu_structure;
                model.max_degree = max(sum(model.structure));
                model.median_degree = median(sum(model.structure));
                model.mean_degree = mean(sum(model.structure));
                model.rms_degree = rms(sum(model.structure));

                self.models{end+1} = model;
            end
        end

        % naive bayes structure
        % structure with no edges
        function self = do_naive_bayes_structure_learning(self)
            fprintf('Naive Bayes structure, i.e. structure with no edges');

            naive_bayes_structure = zeros(size(self.x_train,2));

            for k = 1:numel(self.p_lambda_sequence)
                model = struct();
                model.p_lambda = self.p_lambda_sequence(k);
                model.s_lambda = 0;
                model.density = 0;
                model.structure = naive_bayes_structure;
                model.max_degree = max(sum(model.structure));
                model.median_degree = median(sum(model.structure));
                model.mean_degree = mean(sum(model.structure));
                model.rms_degree = rms(sum(model.structure));

                self.models{end+1} = model;
            end

        end

        % runs the parameter estimation code using the structures learned
        function self = do_parameter_estimation(self, varargin)
            parser = inputParser;
            parser.addParamValue('BCFW_max_iterations', 75000, @isnumeric);
            parser.addParamValue('BCFW_fval_epsilon', .1, @isnumeric);
            parser.addParamValue('compute_true_logZ', true, @islogical);
            parser.addParamValue('reweight_denominator', 'max_degree');

            if islogical(self.x_train) && islogical(self.x_test)
                hidden_model = false;
            else
                hidden_model = true;
            end
            parser.addParamValue('hidden_model', hidden_model, @islogical);

            parser.parse(varargin{:});
            BCFW_max_iterations = parser.Results.BCFW_max_iterations;
            BCFW_fval_epsilon = parser.Results.BCFW_fval_epsilon;
            self.computed_true_logZ = parser.Results.compute_true_logZ;
            self.hidden_model = parser.Results.hidden_model;
            reweight_denominator = parser.Results.reweight_denominator;

            % for each model, lets estimate its parameters
            for i = 1:numel(self.models)
                model = self.models{i};

                fprintf('\nParameter Estimation: s_lambda=%e; density=%e; p_lambda = %e\n',...
                    model.s_lambda, model.density, model.p_lambda);

                % Converts structure + training samples in the overcomplete parametrization
                % necessary to run the parameter estimation code.
                % Only recompute overcomplete struct if structure has
                % changed
                if i == 1 || any(any(self.models{i-1}.structure ~= model.structure))
                    overcomplete_struct = samples_to_overcomplete(self.x_train, model.structure);
                end

                % define reweight parameter
                if ischar(reweight_denominator)
                    if strcmp(reweight_denominator, 'max_degree')
                        reweight = 2/model.max_degree;
                    elseif strcmp(reweight_denominator, 'median_degree')
                        reweight = 2/model.median_degree;
                    elseif strcmp(reweight_denominator, 'mean_degree')
                        reweight = 2/model.mean_degree;
                    elseif strcmp(reweight_denominator, 'rms_degree')
                        reweight = 2/model.rms_degree;
                    else
                        error('Unknown reweighting denominator ''%s''', reweight_denominator);
                    end
                else
                   reweight = 2/reweight_denominator;
                end

                % sanity check
                if reweight > 1
                    reweight = 1;
                end

                % saves the reweight for future cross-validation
                model.reweight = reweight;

                % Create object to do parameter estimation
                loopy_model_train_object = Ising( ...
                    overcomplete_struct.YN, ...
                    overcomplete_struct.YE, ...
                    overcomplete_struct.Ut, ...
                    overcomplete_struct.Vt, ...
                    overcomplete_struct.Ns, ...
                    overcomplete_struct.edges, ...
                    model.p_lambda, ...
                    'checkStuck', false, ...
                    'reweight', reweight);

                % Run Batch C Frank-Wolfe
                bcfw = BCFW(loopy_model_train_object, ...
                    'printInterval', 1000, ...
                    'printTest', 100, ...
                    'printComputedDualityGap', false, ...
                    'MaxIter', BCFW_max_iterations, ...
                    'fvalEpsilon', BCFW_fval_epsilon);
                bcfw.run();

                % get F and G parameters
                model.theta = bcfw.obj.computeParams();

                % Modify F and G


                % compute approximate logZ
                logZ = bcfw.obj.partition_function(model.theta);

                % compute node, edge potentials and adjusted partition
                % function
                fprintf('Converting F and G to node and edge potentials\n');
                [node_pot, edge_pot, logZ_pot] = get_node_and_edge_potentials(model.theta.F,...
                    model.theta.G, logZ, overcomplete_struct.edges{1}');
                model.theta.node_potentials = node_pot;
                model.theta.edge_potentials = edge_pot;
                model.theta.logZ = logZ_pot;

                % compute exact true_logZ and true_node_marginals
                if self.computed_true_logZ
                    fprintf('Starting to run JTA to compute true partition function\n');
                    [true_node_marginals,~,~,~,~,~,true_logZ] = run_junction_tree(node_pot, edge_pot, 'verbose', true);
                    model.theta.true_logZ = true_logZ;
                    model.true_node_marginals = true_node_marginals;
                end

                % Compute likelihoods
                if ~self.hidden_model
                    % Compute training likelihood
                    fprintf('Computing training likelihood\n');
                    model.train_likelihood = compute_avg_log_likelihood( ...
                        model.theta.node_potentials, ...
                        model.theta.edge_potentials, ...
                        model.theta.logZ, ...
                        self.x_train);

                    % Compute test likelihood
                    fprintf('Computing test likelihood\n');
                    model.test_likelihood = compute_avg_log_likelihood( ...
                        model.theta.node_potentials, ...
                        model.theta.edge_potentials, ...
                        model.theta.logZ, ...
                        self.x_test);
                end
                % Compute true test and training likelihood
                if self.computed_true_logZ && self.hidden_model
                    fprintf('Computing true test likelihood\n');
                    model.true_test_likelihood = compute_avg_log_likelihood_hidden( ...
                        model.theta.node_potentials, ...
                        model.theta.edge_potentials, ...
                        model.theta.true_logZ, ...
                        self.x_test);
                elseif self.computed_true_logZ && ~self.hidden_model
                    fprintf('Computing true test likelihood\n');
                    model.true_test_likelihood = compute_avg_log_likelihood( ...
                        model.theta.node_potentials, ...
                        model.theta.edge_potentials, ...
                        model.theta.true_logZ, ...
                        self.x_test);

                    fprintf('Computing true training likelihood\n');
                    model.true_train_likelihood = compute_avg_log_likelihood( ...
                        model.theta.node_potentials, ...
                        model.theta.edge_potentials, ...
                        model.theta.true_logZ, ...
                        self.x_train);
                end

                model.pending_parameter_estimation = false;
                self.models{i} = model;
            end
            fprintf('Finished estimating parameters.\n');
        end


        % return the model (a SingleLoopyModel object) that has the highest
        % test likelihood
        function [best_model] = get_best_model(self)
            % if true likelihood available, use it
            if self.computed_true_logZ
                [~, best_model_index] = max(cellfun(@(m) m.true_test_likelihood, self.models));
            else
                % find the highest interior maxima for each density, and
                % then pick the best
                model_structs = cellfun(@(x) struct(...
                    's_lambda',find(self.s_lambda_sequence==x.s_lambda),...
                    'p_lambda',find(self.p_lambda_sequence==x.p_lambda),...
                    'test_likelihood',x.test_likelihood,...
                    'density',find(self.density_sequence==x.density)...
                    ), self.models);

                best_by_density = struct([]);
                for d = 1:numel(self.density_sequence)
                    density_structs = model_structs([model_structs.density] == d);
                    dimensions_likelihood_matrix = cell2mat(reshape(struct2cell(density_structs), 4,[])');
                    likelihood_grid = full(sparse(dimensions_likelihood_matrix(:,1),...
                        dimensions_likelihood_matrix(:,2),...
                        dimensions_likelihood_matrix(:,3)));

                    best_row = 0;
                    best_col = 0;
                    best_like = -inf;
                    for i = 2:(numel(self.s_lambda_sequence)-1)
                        for j = 2:(numel(self.p_lambda_sequence)-1)
                            if likelihood_grid(i,j) > likelihood_grid(i-1,j) && ...
                                    likelihood_grid(i,j) > likelihood_grid(i,j-1) && ...
                                    likelihood_grid(i,j) > likelihood_grid(i+1,j) && ...
                                    likelihood_grid(i,j) > likelihood_grid(i,j+1) && ...
                                    likelihood_grid(i,j) > likelihood_grid(i-1,j-1) && ...
                                    likelihood_grid(i,j) > likelihood_grid(i-1,j+1) && ...
                                    likelihood_grid(i,j) > likelihood_grid(i+1,j-1) && ...
                                    likelihood_grid(i,j) > likelihood_grid(i+1,j+1)
                                if likelihood_grid(i,j) > best_like
                                    best_row = i;
                                    best_col = j;
                                    best_like = likelihood_grid(i,j);
                                end
                            end
                        end
                    end
                    if best_row ~= 0 && best_col ~= 0
                        best_by_density(end+1).s_lambda = self.s_lambda_sequence(best_row);
                        best_by_density(end).p_lambda = self.p_lambda_sequence(best_col);
                        best_by_density(end).like = best_like;
                        best_by_density(end).density = self.density_sequence(d);
                    end
                end

                if length(best_by_density) > 0
                    [~,best_index] = max([best_by_density.like]);
                    best_s_lambda = best_by_density(best_index).s_lambda;
                    best_p_lambda = best_by_density(best_index).p_lambda;
                    best_density = best_by_density(best_index).density;

                    [~, best_model_index] = max(cellfun(@(m) (m.s_lambda==best_s_lambda) && (m.p_lambda==best_p_lambda) && (m.density==best_density) , self.models));
                else
                   % error('No best model found');
	           fprintf('No best model found\n');
                   best_model_index = 106;
                end
            end

            best_model = SingleLoopyModel( ...
                self.x_train, self.x_test, ...
                self.models{best_model_index}, ...
                self.variable_names);
        end


        %##################################################################
        % Visualization Functions
        %##################################################################
        function visualize_best_model_graph(self)
            best_model_object = self.get_best_model();
            best_model_object.plot_graph();
        end

        function parameter_plots(self)
            best_model_object = self.get_best_model();

            % PLOT 1
            % surf plot showing true test likelihood for each combination
            % of s_lambda and p_lambda, using density of the best model
            likelihoods = zeros(length(self.p_lambda_sequence), length(self.s_lambda_sequence));
            for i = 1:length(self.models)
                model = self.models{i};
                if model.density == best_model_object.density
                    p_lambda_index = find(self.p_lambda_sequence == model.p_lambda);
                    s_lambda_index = find(self.s_lambda_sequence == model.s_lambda);
                    likelihoods(p_lambda_index, s_lambda_index) = model.true_test_likelihood;
                end
            end
            figure(1);
            surf(self.s_lambda_sequence, self.p_lambda_sequence, likelihoods);
            xlabel('Structure learning \lambda');
            ylabel('Parameter estimation \lambda');
            zlabel('True avg. test log-likelihood');

            % PLOT 2
            % p_lambda vs test log likelihood
            true_likelihoods = zeros(1,length(self.p_lambda_sequence));
            approx_likelihoods = zeros(1,length(self.p_lambda_sequence));
            for i = 1:length(self.models)
                model = self.models{i};
                if model.s_lambda == best_model_object.s_lambda && model.density == best_model_object.density
                    p_lambda_index = find(self.p_lambda_sequence == model.p_lambda);
                    true_likelihoods(p_lambda_index) = model.true_test_likelihood;
                    approx_likelihoods(p_lambda_index) = model.test_likelihood;
                end
            end
            figure(2);
            curve1 = semilogx(self.p_lambda_sequence, true_likelihoods);
            set(curve1, 'LineWidth',2); set(curve1, 'Color','b');
            hold on;
            curve2 = semilogx(self.p_lambda_sequence, approx_likelihoods);
            set(curve2, 'LineWidth',2); set(curve2, 'Color','r');
            ylabel('Avg. test log-likelihood');
            xlabel('Parameter estimation \lambda');
            grid on;
            legend('True', 'Approx.');
            hold off;

            % PLOT 3
            % s_lambda vs test log likelihood
            true_likelihoods = zeros(1,length(self.s_lambda_sequence));
            approx_likelihoods = zeros(1,length(self.s_lambda_sequence));
            for i = 1:length(self.models)
                model = self.models{i};
                if model.p_lambda == best_model_object.p_lambda && model.density == best_model_object.density
                    s_lambda_index = find(self.s_lambda_sequence == model.s_lambda);
                    true_likelihoods(s_lambda_index) = model.true_test_likelihood;
                    approx_likelihoods(s_lambda_index) = model.test_likelihood;
                end
            end
            figure(3);
            curve1 = semilogx(self.s_lambda_sequence, true_likelihoods);
            set(curve1, 'LineWidth',2); set(curve1, 'Color','b');
            hold on;
            curve2 = semilogx(self.s_lambda_sequence, approx_likelihoods);
            set(curve2, 'LineWidth',2); set(curve2, 'Color','r');
            ylabel('Avg. test log-likelihood');
            xlabel('Structure learning \lambda');
            grid on;
            legend('True', 'Approx.');
            hold off;

            % PLOT 4
            % density vs test log likelihood
            true_likelihoods = zeros(1,length(self.density_sequence));
            approx_likelihoods = zeros(1,length(self.density_sequence));
            for i = 1:length(self.models)
                model = self.models{i};
                if model.p_lambda == best_model_object.p_lambda && model.s_lambda == best_model_object.s_lambda
                    density_index = find(self.density_sequence == model.density);
                    true_likelihoods(density_index) = model.true_test_likelihood;
                    approx_likelihoods(density_index) = model.test_likelihood;
                end
            end
            figure(4);
            curve1 = semilogx(self.density_sequence, true_likelihoods);
            set(curve1, 'LineWidth',2); set(curve1, 'Color','b');
            hold on;
            curve2 = semilogx(self.density_sequence, approx_likelihoods);
            set(curve2, 'LineWidth',2); set(curve2, 'Color','r');
            ylabel('Avg. test log-likelihood');
            xlabel('Density');
            grid on;
            legend('True', 'Approx.');
            hold off;
        end

        % Plots the approximate likelihood grid for each density, by
        % regularization parameter combinations.
        function plot_likelihood_grid(self)
            % find the highest interior maxima for each density, and
            % then pick the best
            model_structs = cellfun(@(x) struct(...
                's_lambda',find(self.s_lambda_sequence==x.s_lambda),...
                'p_lambda',find(self.p_lambda_sequence==x.p_lambda),...
                'test_likelihood',x.test_likelihood,...
                'density',find(self.density_sequence==x.density)...
                ), self.models);

            for d = 1:numel(self.density_sequence)
                density_structs = model_structs([model_structs.density] == d);
                dimensions_likelihood_matrix = cell2mat(reshape(struct2cell(density_structs), 4,[])');
                likelihood_grid = full(sparse(dimensions_likelihood_matrix(:,1),...
                    dimensions_likelihood_matrix(:,2),...
                    dimensions_likelihood_matrix(:,3)));

                figure(40+floor(self.density_sequence(d)*100));
                colormap bone;
                contourf(self.p_lambda_sequence, self.s_lambda_sequence, likelihood_grid, 40);
                xlabel('\lambda_p'); ylabel('\lambda_s');
                set(gca,'xscale','log'); set(gca,'yscale','log');
                set(gca, 'XTick', self.p_lambda_sequence); set(gca, 'YTick', self.s_lambda_sequence);
                set(gca, 'XTickLabel', sprintf('%.0e| |',self.p_lambda_sequence(1:2:end))); set(gca, 'YTickLabel', sprintf('%.3f|',self.s_lambda_sequence));
                set(gca,'XMinorTick','off','YMinorTick','off');
                title(sprintf('True Log-likelihood (density: %.2f)', self.density_sequence(d)));
            end
        end

        % Merge Models
        % This method exists because instead of training all models in a
        % single machine, we split the models in multiple
        % LoopyModelCollection objects across many machines. This method
        % merges two LoopyModelCollection objects into one, keeping the
        % models ordered by: s_lambda, density, p_lambda
        function self = merge_model_collections(self, other_collection)
            % concatenate cellarray
            merged_models = [self.models other_collection.models];

            % create struct array with fields s_lambda, density and
            % p_lambda. Then convert struct array into a matrix. We end up
            % with a matrix with one model per row, and three columns:
            % s_lambda, density and p_lambda
            models_params = cellfun(@(x) struct('s_lambda',x.s_lambda, 'density',x.density, 'p_lambda',x.p_lambda), merged_models);
            models_params = reshape(struct2cell(models_params), 3, [])';

            % now sort model_params by columns 1,2 and 3 respectively, and
            % get back a index vector
            [~, sorting_indexes] = sortrows(models_params, [1 2]);

            % sort the models
            self.models = merged_models(sorting_indexes);

            % update density_sequence, s_lambda_sequence, p_lambda_sequence
            densities = cellfun(@(x) x.density, self.models);
            densities = sort(unique(densities));
            self.density_sequence = densities;
            s_lambdas = cellfun(@(x) x.s_lambda, self.models);
            s_lambdas = sort(unique(s_lambdas));
            self.s_lambda_sequence = s_lambdas;
            p_lambdas = cellfun(@(x) x.p_lambda, self.models);
            p_lambdas = sort(unique(p_lambdas));
            self.p_lambda_sequence = p_lambdas;
        end

        function self = set_p_lambdas(self, p_lambda_count, p_lambda_min, p_lambda_max)
            p_lambda_min_exp = log10(p_lambda_min);
            p_lambda_max_exp = log10(p_lambda_max);
            self.p_lambda_sequence = logspace(p_lambda_min_exp, p_lambda_max_exp, p_lambda_count);

            for i = 1:numel(self.s_lambda_sequence)
                for j = 1:numel(self.density_sequence)
                    for k = 1:numel(self.p_lambda_sequence)
                        index = (i-1)*length(self.density_sequence)*length(self.p_lambda_sequence)...
                            + (j-1)*length(self.p_lambda_sequence) + k;
                        model = self.models{index};
                        model.p_lambda = self.p_lambda_sequence(k);
                        self.models{index} = model;
                    end
                end
            end
        end
    end

    methods(Static)
        function parser = parse_and_validate_input(x_train, x_test, varargin)
            parser = inputParser;
%             is_logical_matrix = @(x) ismatrix(x) && islogical(x) && ~isempty(x);
            is_logical_matrix = @(x) ismatrix(x) && ~isempty(x);
            parser.addRequired('x_train', is_logical_matrix);
            parser.addRequired('x_test', is_logical_matrix);

            is_positive = @(x) x >= 0;
            parser.addParamValue('s_lambda_count', 10, is_positive);
            parser.addParamValue('s_lambda_min', 1e-04, is_positive);
            parser.addParamValue('s_lambda_max', 0.9, is_positive);
            parser.addParamValue('density_count', 2, is_positive);
            parser.addParamValue('density_min', 0.05, @isscalar);
            parser.addParamValue('density_max', 0.09, @isscalar);
            parser.addParamValue('p_lambda_count', 20, is_positive);
            parser.addParamValue('p_lambda_min', 1e-03, is_positive);
            parser.addParamValue('p_lambda_max', 1e09, is_positive);

            % default value for variable names are just sequential numbers
            sequential_number_strings = regexp(num2str(1:(size(x_train,2))),'\s+','split');
            parser.addParamValue('variable_names', sequential_number_strings, @iscellstr);

            parser.parse(x_train, x_test, varargin{:});
        end
    end


end
