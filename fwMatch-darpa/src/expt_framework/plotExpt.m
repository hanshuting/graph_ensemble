function plotExpt(name, id)
% plotExpt  Run the plot.m file in expt, if it exists.

    exptDir = sprintf('expt/%s', name);
    
    addpath(exptDir);
    
    % Ingest the params struct by evaluating the config file.
    config = sprintf('config%d', id);
    eval(config);
    
    % Automatically log everything that comes onscreen henceforth.
    %diary(sprintf('%s/diary%d.txt', exptDir, id));
    
    % Inject variables into params.
    params.saveFile       = sprintf('%s/result%d.mat', exptDir, id);
    params.checkpointFile = sprintf('%s/%s_cp.mat',      exptDir, config);
    params.logFile        = sprintf('%s/%s_log.txt',     exptDir, config);

    params.figPrefix      = sprintf('%s/result%d_fig_', exptDir, id);

    results = load(sprintf('%s/result%d.mat', exptDir, id));
    
    % Run!
    eval('plot_results(results, params)');    
    
    % Clean up path, in case we were called interactively.
    rmpath(exptDir);
end

