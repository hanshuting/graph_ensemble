function create_config_files(varargin)
    % Remove old config files
    system('rm -f config*.m');
    system('rm -f yeti_config.sh');

    % guarantees folders exist
    system('mkdir -p job_logs');
    system('mkdir -p results');
    system('mkdir -p yeti_logs');

    parser = inputParser;

    parser.addParameter('datapath', [], @ischar);
    parser.addParameter('experiment_name', 'spikes_opto_on', @ischar);
    parser.addParameter('email_for_notifications', 'jds2270@columbia.edu', @ischar);
    parser.addParameter('yeti_user', 'jds2270', @ischar);

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

    parser.addParameter('time_span', 1, @isscalar);

    parser.addParameter('num_shuffle', 100, @isscalar);

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

    time_span = parser.Results.time_span;

    num_shuffle = parser.Results.num_shuffle;

     display('CHECK INFO BELOW');
    display(sprintf('Writing config files for yeti user %s', yeti_user));
    display(sprintf('Experiment name: %s', experiment_name));
    display(sprintf('Total jobs to be submitted: %d', num_shuffle));

    % Write config files
    for config_file_count=1:num_shuffle
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
            fprintf(fid,'[params.data, params.variable_names, params.stimuli] = get_real_data(%d);\n',config_file_count);
        else
            fprintf(fid,'[params.data, params.variable_names, params.stimuli] = get_dataset(''%s%d%s'');\n', ...
                    datapath(1:end-4), int2str(config_file_count), datapath(end-3:end));
        end

        if strcmp(structure_type, 'loopy')
            % slambda
            fprintf(fid,'s_lambdas = logspace(%f,%f,%d);\n', log10(s_lambda_min), log10(s_lambda_max), s_lambdas_per_split*s_lambda_splits);
            fprintf(fid,'params.s_lambda_count = %d;\n', s_lambdas_per_split);
            fprintf(fid,'params.s_lambda_min = s_lambdas(%d);\n', 1);
            fprintf(fid,'params.s_lambda_max = s_lambdas(%d);\n', s_lambdas_per_split);
            % density
            fprintf(fid,'densities = linspace(%f,%f,%d);\n', density_min, density_max, densities_per_split*density_splits);
            fprintf(fid,'params.density_count = %d;\n', densities_per_split);
            fprintf(fid,'params.density_min = densities(%d);\n', 1);
            fprintf(fid,'params.density_max = densities(%d);\n', densities_per_split);

            % lookback time span
            fprintf(fid,'params.time_span = %d;\n', time_span);
        end

        % plambda
        fprintf(fid,'p_lambdas = logspace(%f,%f,%d);\n', log10(p_lambda_min), log10(p_lambda_max), p_lambdas_per_split*p_lambda_splits);
        fprintf(fid,'params.p_lambda_count = %d;\n', p_lambdas_per_split);
        fprintf(fid,'params.p_lambda_min = p_lambdas(%d);\n', 1);
        fprintf(fid,'params.p_lambda_max = p_lambdas(%d);\n', p_lambdas_per_split);

        fclose(fid);
    end

    % Write YETI script
    fid = fopen('yeti_config.sh','w');

    fprintf(fid,'#!/bin/sh\n');
    fprintf(fid,'#yeti_config.sh\n\n');
    fprintf(fid,'#Torque script to run Matlab program\n');

    fprintf(fid,'\n#Torque directives\n');
    fprintf(fid,'#PBS -N %s\n', experiment_name);
    fprintf(fid,'#PBS -W group_list=yetibrain\n');
    if(strcmp( structure_type, 'loopy') == 1)
      fprintf(fid,'#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=8000mb\n');
    else
      fprintf(fid,'#PBS -l nodes=1:ppn=1,walltime=2:00:00,mem=8000mb\n');
    end
    fprintf(fid,'#PBS -m af\n');
    fprintf(fid,'#PBS -M %s\n', email_for_notifications);
    fprintf(fid,'#PBS -V\n');
    fprintf(fid,'#PBS -t 1-%d\n',num_shuffle);

    expt_dir = pwd;
    run_dir = expt_dir(1:regexp(expt_dir, 'fwMatch-darpa/','end'));
    fprintf(fid,'\n#set output and error directories (SSCC example here)\n');
    fprintf(fid,'#PBS -o localhost:%s/yeti_logs/\n', expt_dir);
    fprintf(fid,'#PBS -e localhost:%s/yeti_logs/\n', expt_dir);

    fprintf(fid,'\n#Command below is to execute Matlab code for Job Array (Example 4) so that each part writes own output\n');
    fprintf(fid,'cd %s\n', run_dir);
    fprintf(fid,'./run.sh %s $PBS_ARRAYID > expt/%s/job_logs/matoutfile.$PBS_ARRAYID\n', experiment_name, experiment_name);
    fprintf(fid,'#End of script\n');

fclose(fid);

% Kate: we should write the start_jobs.sh script, too, since it includes the experiment name var
    fid = fopen('start_jobs.sh','w');

    fprintf(fid,'#If you are having problem with line endings use ":set ff=unix" in vim\n\n');

%Kate - I do not want to remove any of the already computed model files; do that manually
fprintf(fid,'rm -f ./results/result*.mat\n');
%fprintf(fid,'rm -f ./results/model_calcium.mat\n');
fprintf(fid,'rm -f ./yeti_logs/*\n');
fprintf(fid,'rm -f ./job_logs/*\n');

fprintf(fid,'cd ../.. && qsub expt/%s/yeti_config.sh\n', experiment_name);
fclose(fid);
%Kate: make sure they are executable:
system('chmod +x start_jobs.sh yeti_config.sh');

end

