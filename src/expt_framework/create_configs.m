function create_configs(varargin)
% Expects to be run in the working directory
    % Remove old config files
    system('rm -f config*.m');
    system('rm -f yeti_config.sh');

    % guarantees folders exist
    system('mkdir -p job_logs');
    system('mkdir -p results');
    system('mkdir -p yeti_logs');

    parser = inputParser;

    parser.addParameter('datapath', [], @ischar);
    parser.addParameter('experiment_name', 'experiment', @ischar);
    parser.addParameter('email_for_notifications', 'UNI@columbia.edu', @ischar);
    parser.addParameter('yeti_user', 'UNI', @ischar);

    parser.addParameter('training_test_split', .8, @isscalar);
    parser.addParameter('BCFW_max_iterations', 75000, @isscalar);
    parser.addParameter('structure_type', 'loopy', @ischar);
    parser.addParameter('compute_true_logZ', false, @islogical);
    parser.addParameter('reweight_denominator', 'max_degree');

    parser.addParameter('s_lambda_splits', 1, @isscalar);
    parser.addParameter('s_lambdas_per_split', 1, @isscalar);
    parser.addParameter('s_lambda_min', 1e-01, @isscalar);
    parser.addParameter('s_lambda_max', 1e-01, @isscalar);

    parser.addParameter('density_splits', 1, @isscalar);
    parser.addParameter('densities_per_split', 1, @isscalar);
    parser.addParameter('density_min', 0.05, @isscalar);
    parser.addParameter('density_max', 0.05, @isscalar);

    parser.addParameter('p_lambda_splits', 1, @isscalar);
    parser.addParameter('p_lambdas_per_split', 1, @isscalar);
    parser.addParameter('p_lambda_min', 1e+01, @isscalar);
    parser.addParameter('p_lambda_max', 1e+01, @isscalar);

    parser.addParameter('edges', 'simple', @ischar);
    parser.addParameter('no_same_neuron_edges', true, @islogical);
    parser.addParameter('time_span', 1, @isscalar);

    parser.parse(varargin{:})

    datapath = parser.Results.datapath;
    experiment_name = parser.Results.experiment_name;
    email_for_notifications = parser.Results.email_for_notifications;
    yeti_user = parser.Results.yeti_user;

    training_test_split = parser.Results.training_test_split;
    BCFW_max_iterations = parser.Results.BCFW_max_iterations;
    structure_type = parser.Results.structure_type;
    compute_true_logZ = parser.Results.compute_true_logZ;
    if (compute_true_logZ); compute_true_logZ_str='true'; else;  compute_true_logZ_str='false'; end
    reweight_denominator = parser.Results.reweight_denominator;

    s_lambda_splits = parser.Results.s_lambda_splits;
    s_lambdas_per_split = parser.Results.s_lambdas_per_split;
    s_lambda_min = parser.Results.s_lambda_min;
    s_lambda_max = parser.Results.s_lambda_max;

    density_splits = parser.Results.density_splits;
    densities_per_split = parser.Results.densities_per_split;
    density_min = parser.Results.density_min;
    density_max = parser.Results.density_max;

    p_lambda_splits = parser.Results.p_lambda_splits;
    p_lambdas_per_split = parser.Results.p_lambdas_per_split;
    p_lambda_min = parser.Results.p_lambda_min;
    p_lambda_max = parser.Results.p_lambda_max;

    edges = parser.Results.edges;
    no_same_neuron_edges = parser.Results.no_same_neuron_edges;
    if (no_same_neuron_edges); no_same_neuron_edges_str='true'; else;  no_same_neuron_edges_str='false'; end
    time_span = parser.Results.time_span;

     display('CHECK INFO BELOW');
    if ~strcmp(yeti_user, "UNI")
        display(sprintf('Writing config files for yeti user %s', yeti_user));
    end
    display(sprintf('Experiment name: %s', experiment_name));
    display(sprintf('Total configs files to be created: %d', p_lambda_splits*s_lambda_splits*density_splits));

    % Write config files
    config_file_count = 0;
    for i=1:p_lambda_splits
        for j=1:s_lambda_splits
            for k=1:density_splits
                config_file_count = config_file_count + 1;
                fid = fopen(sprintf('config%d.m',config_file_count),'w');
                fprintf(fid,'params.split = %f;\n', training_test_split);
                fprintf(fid,'params.BCFW_max_iterations = %d;\n', BCFW_max_iterations);
                fprintf(fid,'params.structure_type = ''%s'';\n', structure_type);
                fprintf(fid,'params.compute_true_logZ = %s;\n', compute_true_logZ_str);
                if ischar(reweight_denominator)
                    fprintf(fid,'params.reweight_denominator = ''%s'';\n', reweight_denominator);
                else
                    fprintf(fid,'params.reweight_denominator = %d;\n', reweight_denominator);
                end

                % get real data (params.data)
                if isempty(datapath)
                    fprintf(fid,'[params.data, params.variable_names, params.stimuli] = get_real_data();\n');
                else
                    fprintf(fid,'[params.data, params.variable_names, params.stimuli] = get_dataset(''%s'');\n', datapath);
                end

                if strcmp(structure_type, 'loopy')
                    % slambda
                    fprintf(fid,'s_lambdas = logspace(%f,%f,%d);\n', log10(s_lambda_min), log10(s_lambda_max), s_lambdas_per_split*s_lambda_splits);
                    fprintf(fid,'params.s_lambda_count = %d;\n', s_lambdas_per_split);
                    fprintf(fid,'params.s_lambda_min = s_lambdas(%d);\n', (j-1)*s_lambdas_per_split + 1);
                    fprintf(fid,'params.s_lambda_max = s_lambdas(%d);\n', j*s_lambdas_per_split);

                    % density
                    fprintf(fid,'densities = linspace(%f,%f,%d);\n', density_min, density_max, densities_per_split*density_splits);
                    fprintf(fid,'params.density_count = %d;\n', densities_per_split);
                    fprintf(fid,'params.density_min = densities(%d);\n', (k-1)*densities_per_split + 1);
                    fprintf(fid,'params.density_max = densities(%d);\n', k*densities_per_split);

                    % edge constraint parameters
                    fprintf(fid,'params.edges = ''%s'';\n', edges);
                    fprintf(fid,'params.no_same_neuron_edges = %s;\n', no_same_neuron_edges_str);

                    % lookback time span
                    fprintf(fid,'params.time_span = %d;\n', time_span);
                end

                % plambda
                fprintf(fid,'p_lambdas = logspace(%f,%f,%d);\n', log10(p_lambda_min), log10(p_lambda_max), p_lambdas_per_split*p_lambda_splits);
                fprintf(fid,'params.p_lambda_count = %d;\n', p_lambdas_per_split);
                fprintf(fid,'params.p_lambda_min = p_lambdas(%d);\n', (i-1)*p_lambdas_per_split + 1);
                fprintf(fid,'params.p_lambda_max = p_lambdas(%d);\n', i*p_lambdas_per_split);

                fclose(fid);
            end
        end
    end

end

